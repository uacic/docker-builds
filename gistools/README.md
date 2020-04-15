# GIS_tools
A collection of python scripts for various GIS tasks. These scripts were developed for data processing for the world's largest phenotyping platform at the University of Arizona's Maricopa Agricultural Center. 

### Docker run
```
docker run --rm -v $(pwd):/mnt acicarizona/gistools --csv "/mnt/2020-01-08.csv" "/mnt/gpscorrect_out"
```
## Editing GPS coordinates of TIF images
The script `edit_gps.py` edits the corner coordinates within your TIF file. Just feed it the directory where your files are located and a CSV file with your coordinates which must contain the headers: Filename, Upper left, Lower right.
