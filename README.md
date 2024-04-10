# WG26 Connectathon Consumer Code

This is the source code for the IDC's entry into the WG26 Connectathon on Bulk
Microscopy Simple Annotations (ANNs).

It is a simple Python script intended to demonstrate how to leverage existing
open source Python libraries (including highdicom, dicomweb-client, shapely,
and numpy) to enable parsing and simple computations on information stored in
annotation objects.

The script pulls all ANN objects from it finds in both archives used for the
connectathon, loops over the annotations present and calculates the following
properties of each annotation:

- Area
- Height
- Width
- Perimeter
- Centroid

It then puts these results into a CSV file along with some other information
about the annotation (its label, number and property type and code). There is
single CSV file output for each copy of an ANN object in a given archive.
Each CSV contains a single line per annotation.

The entire process happens over DICOMweb and there is no need to write out
DICOM objects to the filesystem.

### Installation/Use

Use the provided requirements file to install the required packages into a
virtualenv, then just execute the script with no arguments.
