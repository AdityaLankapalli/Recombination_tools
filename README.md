# Recombination_tools
Scripts to analyse recombination sites identified in bacteria

## CFML_sitesremover.py  
Removes the sites of recombination identified by ClonalFrameML from an alignment file. 
The program removes sites mentioned in the importation_status file from ClonalFrameML.

#### How to run

```
python3 CFML_sitesremover.py -f <Recomb sites> -i <Input fasta> -o <Output fasta>
```

###### Options

-f    Recombination sites specified in importation_status.txt of ClonalFrameML <br>
-i    Input fasta filename <br>
-o    Output fasta filename <br>

