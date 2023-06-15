Calculate the weight
$ ../get_wts.py 6v8x.pdb.gz A CD -o 6v8x.A_CD.txt

Calculate the epitope distance
$ ../ep_dist.py -msa V703.with_ref.fasta -w 6v8x.A_CD.txt -p "V703" -r "ref.C|VRC01" -o V703.epdist.txt
$ ../ep_dist.py -msa V704.with_ref.fasta -w 6v8x.A_CD.txt -p "V704" -r "ref.B|VRC01" -o V704.epdist.txt
