# BioPAN - A web-based tool to explore mammalian lipidome metabolic pathways
BioPAN for Bioinformatics Methodology For Pathway Analysis is a tool that allows users to upload their own mammalian lipidomics dataset and perform a pathway analysis. In this analysis, you can explore systematic changes in lipid pathways at different levels: lipid subclass and lipid molecular species. Those pathways will be highlighted and changes in gene activity will also be predicted.


## Getting Started
Install VirtualBox (https://www.virtualbox.org/) and create a new virtual machine with the available BioPAN.ova file (File > Import Appliance).
Make sure that in "Settings > Network > Advanced > Port forwarding" you have the host port 8080.
You can enter the VM with username: root and password: access_biopan.
BioPAN is available at: http://127.0.0.1:8080/


### Documentation
Documentation on how BioPAN works is available at http://127.0.0.1:8080/doc/
Documentation on BioPAN code is available under "biopan > doc > code". There is documentations on: the database, R scripts and web structure.


### Directory and file structure

```bash
├── biopan
|   ├── css
│   │   ├── *.css
|   ├── database
│   │   ├── *.sql
│   ├── doc
│   │   ├── _images
│   │   ├── _static
│   │   ├── code
│   │   │   ├── db_diagram.pdf
│   │   │   ├── rScripts_descr.pdf
│   │   │   ├── web_structure.pdf
│   │   ├── *.html
│   ├── images
│   ├── js
│   │   ├── librairies
│   │   ├── *.js
│   ├── R
│   │   ├── compute_pathway.r
│   │   ├── filter_pathway.r
│   │   ├── lib_parse_data.r
│   │   ├── lib_pathway_analysis.r
│   │   ├── lib_process_data.r
│   │   ├── lib.r
│   │   ├── parse_data.r
│   │   ├── process_data.r
│   ├── resources
│   │   ├── lipid_nodes.csv
│   │   ├── sample.csv
│   │   ├── smallSample.csv
│   ├── stylesheets
│   │   ├── *.cycss
│   ├── *.php
│   ├── runLipidLynxX.sh

```

