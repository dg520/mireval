from __future__ import division
import shlex, subprocess
import commands

root="/home/shared"

cmd="ls /home/shared/mireval/genomes/"
sta,out=commands.getstatusoutput(cmd)
species=out.split("\n")

for sp in species:

    if sp=="homo_sapiens":                        # GRCh37/hg19
        method_avbl=(1,1,1,1)        
    if sp=="mus_musculus":                        # GRCm38/mm10
        method_avbl=(1,1,1,1)
    if sp=="bos_taurus":                          # Btau_4.6.1/bosTau7
        method_avbl=(0,0,1,1)
    if sp=="caenorhabditis_briggsae":             # cb3
        method_avbl=(0,0,1,1)
    if sp=="caenorhabditis_elegans":              # WS220/ce10
        method_avbl=(1,0,1,1)
    if sp=="canis_familiaris":                    # canFam2; not up-to-date
        method_avbl=(1,0,1,0)
    if sp=="drosophila_melanogaster":             # BDGP_5.0/dm3
        method_avbl=(1,0,1,1)
    if sp=="danio_rerio":                         # Zv9/danRer7
        method_avbl=(1,1,1,0)
    if sp=="gallus_gallus":                       # galGal4/ICGSC Gallus_gallus-4.0
        method_avbl=(0,0,1,1)
    if sp=="xenopus_tropicalis":                  # xenTro3/JGI_4.2
        method_avbl=(1,0,1,0)
    if sp=="rattus_norvegicus":                   # rn4/RGSC_3.4; not up-to-date
        method_avbl=(1,0,1,1)
    if sp=="pan_troglodytes":                     # CSAC_2.1.4/panTro4
        method_avbl=(0,0,1,1)
    if sp=="macaca_mulatta":                      # rheMac2/MMUL_1.0; not up-to-date
        method_avbl=(0,0,1,1)
    if sp=="monodelphis_domestica":               # MonDom5
        method_avbl=(1,0,1,1)
    if sp=="fugu_rubripes":                       # FUGU5/fr3
        method_avbl=(1,0,1,0)

    if sp=="pongo_albelii":                       # ponAbe2/WUSTL_2.0.2
        method_avbl=(1,0,0,0)
    if sp=="gorilla_gorilla":                     # gorGor3.1/gorGor3
        method_avbl=(1,0,1,0)
    if sp=="branchiostoma_floridae":              # braFlo1/JGI_1.0
        method_avbl=(1,0,0,0)
    if sp=="petromyzon_marinus":                  # WUGSC 7.0/petMar2
        method_avbl=(1,0,1,0)
    if sp=="oryzias_latipes":                     # oryLat2/MEDAKA1
        method_avbl=(1,0,1,0)
    if sp=="geospiza_fortis":                     # GeoFor_1.0/geoFor1
        method_avbl=(1,1,0,0)
    if sp=="ornithorhynchus_anatinus":            # ornAna1/WUSTL_5.0.1
        method_avbl=(1,0,1,1)
    if sp=="gasterosteus_aculeatus":              # gasAcu1/BI_1.0
        method_avbl=(1,0,0,0)
    if sp=="saccharomyces_cerevisiae":            # sacCer3/Apr2011
        method_avbl=(1,0,0,1)

    if sp=="sus_scrofa":                          # SGSC_Sscrofa10.2/susScr3
        method_avbl=(0,0,1,1)
    if sp=="felis_catus":                         # Felis_catus_6.2/felCat5
        method_avbl=(0,0,0,0)
    if sp=="oryctolagus_cuniculus":               # Broad/oryCun2
        method_avbl=(0,0,0,1)
    if sp=="ovis_aries":                          # Ovis_aries_1.0/oviAri1
        method_avbl=(0,0,1,1)
    if sp=="equus_caballus":                      # equCab2
        method_avbl=(0,0,1,1)
    if sp=="macropus_eugenii":                    # TWGS_Meug_1.1/macEug2
        method_avbl=(0,0,0,0)
    if sp=="sarcophilus_harrisii":                # Devil_ref v7.0/sarHar1
        method_avbl=(0,0,1,0)



    if method_avbl[2]==1 and sp!="rfam":
        f=open(root+"/mireval/genomes/"+sp+"/mir/mir.gff3","r")
        mir_info=f.readlines()
        f.close()
        mir_length=[]
        for mir in mir_info[:-1]:
            if ("primary" in mir) == True and mir[0]!="#":
                t=mir.split("\t")
                mir_length.append(int(t[4])-int(t[3])+1)
        len_mean=int(sum(mir_length)/len(mir_length))
        print sp, len_mean
