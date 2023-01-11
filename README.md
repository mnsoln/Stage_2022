# Stage_2022

Mes scripts de mon stage volontaire à l'IBMP. / My scripts from my internship at IBMP.

**Theme :** ***Analyse des données de séquençage d'A. thaliana. Mise en place d'un pipeline en python pour la détection des niveaux de méthylation. Analyse graphique des données sur R.***


## Installation

```conda create -n methylation python=3.9```

```pip3 install tqdm pandas pysam biopython```


## Description

Diagramme explicatif : https://drive.google.com/file/d/1euW_eSfDgnDULG62Csubbuf5YOSQ8eeX/view?usp=sharing

Uses argparse for arguments.

*Exemple :*

``` python3.9 pos.py --fasta /exemple/test/fasta.fasta --gff /exemple/test/genes_transposons.gff ```

``` python3.9 data.py --reference /exemple/test/ReferencePosC.txt --bismark /exemple/test/bismark.CX_report.txt --debug ```



## Resultats

Pos.py : txt files

Data.py : txt files

Use the rmarkdown file to get graphs (the functions come from the R file). You can also get tables from the R file.
