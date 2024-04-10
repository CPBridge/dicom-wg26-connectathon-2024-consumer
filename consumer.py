#!/usr/bin/env python

from pathlib import Path
import numpy as np
import shapely
import pandas as pd
import pydicom
import highdicom as hd
from scipy.special import ellipe
from dicomweb_client.api import DICOMwebClient
from shapely.measurement import length


if __name__ == "__main__":
    output_dir = Path("./results")
    output_dir.mkdir(exist_ok=True)

    clients = {
        "google": DICOMwebClient(
            url="https://dicomwebproxy-bqmq3usc3a-uc.a.run.app/dicomWeb",
        ),
        "j4care": DICOMwebClient(
            url="https://development.j4care.com:11443/dcm4chee-arc/aets/DCM4CHEE/rs",
        )
    }

    # Loop over clients
    for client_name, client in clients.items():
        print(client_name)
        studies = [pydicom.Dataset.from_json(s) for s in client.search_for_studies()]
        client_dir = output_dir / client_name
        client_dir.mkdir(exist_ok=True)

        # Loop over studies
        for study in studies:

            # Get a list of all annotations in this study
            ann_series = [
                pydicom.Dataset.from_json(s)
                for s in client.search_for_series(
                    study_instance_uid=study.StudyInstanceUID,
                    search_filters={"Modality": "ANN"}
                )
            ]

            for s in ann_series:
                ann_dcms = client.retrieve_series(
                    study_instance_uid=study.StudyInstanceUID,
                    series_instance_uid=s.SeriesInstanceUID,
                )

                for dcm in ann_dcms:
                    ann = hd.ann.MicroscopyBulkSimpleAnnotations.from_dataset(dcm)

                    if ann.annotation_coordinate_type == hd.ann.AnnotationCoordinateTypeValues.SCOORD3D:
                        # There are some 3D annotations but they are out of scope
                        continue

                    # Find and pull metadata for the associated image. This is needed for spacing information
                    ref_ins_uid = ann.ReferencedSeriesSequence[0].ReferencedInstanceSequence[0].ReferencedSOPInstanceUID
                    ref_series_uid = ann.ReferencedSeriesSequence[0].SeriesInstanceUID

                    try:
                        im_meta_json = client.retrieve_instance_metadata(
                            study_instance_uid=study.StudyInstanceUID,
                            series_instance_uid=ref_series_uid,
                            sop_instance_uid=ref_ins_uid,
                        )
                    except:
                        print("Failed to retrieve referenced image.")
                        continue
                    im_meta = pydicom.Dataset.from_json(im_meta_json)

                    spacing_xy = (
                        im_meta.SharedFunctionalGroupsSequence[0]
                        .PixelMeasuresSequence[0]
                        .PixelSpacing
                    )
                    spacing_y, spacing_x = spacing_xy
                    pixel_area = spacing_y * spacing_x
                    assert abs(spacing_x - spacing_y) < 1e-4
                    spacing = spacing_x

                    all_objects = []

                    for grp in ann.get_annotation_groups():
                        graphic_data = grp.get_graphic_data("2D")

                        # Loop over individual annotations in this group as
                        # numpy arrays of (x, y) coordinates
                        for i, coords in enumerate(graphic_data):

                            object_results = {
                                "graphic_type": grp.graphic_type.value,
                                "annotation_group_label": grp.label,
                                "annotation_group_number": grp.number,
                                "annotation_group_property_type_scheme": grp.annotated_property_type.scheme_designator,
                                "annotation_group_property_type_code": grp.annotated_property_type.value,
                                "annotation_group_property_type_meaning": grp.annotated_property_type.meaning,
                                "annotation_group_property_category_scheme": grp.annotated_property_category.scheme_designator,
                                "annotation_group_property_category_code": grp.annotated_property_category.value,
                                "annotation_group_property_category_meaning": grp.annotated_property_category.meaning,
                                "annotation_number": i,
                            }

                            if grp.graphic_type == hd.ann.GraphicTypeValues.POINT:

                                object_results["area_mm2"] = 0.0
                                object_results["height_mm"] = 0.0
                                object_results["width_mm"] = 0.0
                                object_results["boundary_length_mm"] = 0.0
                                object_results["centroid_x_pix"] = coords[0, 0]
                                object_results["centroid_y_pix"] = coords[0, 1]

                            elif grp.graphic_type in (
                                hd.ann.GraphicTypeValues.POLYGON,
                                hd.ann.GraphicTypeValues.RECTANGLE,
                            ):
                                polygon = shapely.Polygon(coords)

                                minx, miny, maxx, maxy = polygon.bounds

                                object_results["area_mm2"] = polygon.area * pixel_area
                                object_results["height_mm"] = (maxy - miny) * spacing
                                object_results["width_mm"] = (maxx - minx) * spacing
                                object_results["boundary_length_mm"] = polygon.length * spacing
                                object_results["centroid_x_pix"] = polygon.centroid.x
                                object_results["centroid_y_pix"] = polygon.centroid.y

                            elif grp.graphic_type == hd.ann.GraphicTypeValues.ELLIPSE:

                                maxx = coords[:, 0].max()
                                minx = coords[:, 0].min()
                                maxy = coords[:, 1].max()
                                miny = coords[:, 1].min()

                                major_axis_length = np.sqrt(
                                    ((coords[0, :] - coords[1, :]) ** 2).sum()
                                ) * spacing
                                minor_axis_length = np.sqrt(
                                    ((coords[2, :] - coords[3, :]) ** 2).sum()
                                ) * spacing

                                eccentricity_sq = 1 - (minor_axis_length ** 2) / (major_axis_length ** 2)

                                object_results["area_mm2"] = np.pi * major_axis_length * minor_axis_length

                                # This is slightly off but not meaningfully so
                                object_results["height_mm"] = (maxy - miny) * spacing
                                object_results["width_mm"] = (maxx - minx) * spacing

                                object_results["boundary_length_mm"] = 4 * major_axis_length * ellipe(eccentricity_sq)

                                centroid = (coords[0, :] + coords[1, :]) / 2.0
                                object_results["centroid_x_pix"] = centroid[0].item()
                                object_results["centroid_y_pix"] = centroid[1].item()


                            elif grp.graphic_type == hd.ann.GraphicTypeValues.POLYLINE:

                                maxx = coords[:, 0].max()
                                minx = coords[:, 0].min()
                                maxy = coords[:, 1].max()
                                miny = coords[:, 1].min()

                                diff = np.diff(coords, axis=0)
                                distances = np.sqrt(diff ** 2).sum(axis=1)
                                length = distances.sum()

                                object_results["area_mm"] = np.nan
                                object_results["height_mm"] = (maxy - miny) * spacing
                                object_results["width_mm"] = (maxx - minx) * spacing
                                object_results["boundary_length_mm"] = length * spacing
                                object_results["centroid_x_pix"] = coords[:, 0].mean()
                                object_results["centroid_y_pix"] = coords[:, 1].mean()

                            all_objects.append(object_results)


                    results_df = pd.DataFrame(all_objects)
                    out_file_name = client_dir / f"{ann.Manufacturer}_{ann.PatientID}_{ann.StudyInstanceUID}_{ann.SeriesInstanceUID}_results.csv"
                    results_df.to_csv(out_file_name)
                    print("Written", str(out_file_name))
