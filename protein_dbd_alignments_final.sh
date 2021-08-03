#AP1
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.2.2_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.1.2_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.4_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.2.1_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.1.1_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.3.2_mammalia_dbd_fasta.fasta
#missing DBD
echo '>Mus_musculus_MAFG-annot_ma-DBD_ma' >> 1.1.3.2_mammalia_dbd_fasta.fasta
echo 'QLKQRRRTLKNRGYAASCRVKRVTQKEELEKQKAELQQEVEKLASENASMKLELDALRSKYEAL' >> 1.1.3.2_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.3.1_mammalia_dbd_fasta.fasta

#twist
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.2.3.2_mammalia_dbd_fasta.fasta
echo '>Mus_musculus_TWIST2-annot_ma-DBD_ma' >> 1.2.3.2_mammalia_dbd_fasta.fasta
echo 'QRILANVRERQRTQSLNEAFAALRKIIPTLPSDKLSKIQTLKLAARYIDFLY' >> 1.2.3.2_mammalia_dbd_fasta.fasta
echo '>Mus_musculus_HAND2-annot_ma-DBD_ma' >> 1.2.3.2_mammalia_dbd_fasta.fasta
echo 'RRGTANRKERRRTQSINSAFAELRECIPNVPADTKLSKIKTLRLATSYIAYLM' >> 1.2.3.2_mammalia_dbd_fasta.fasta
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.3.16_mammalia_dbd_fasta.fasta

#gr
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.1.1.1_mammalia_dbd_fasta.fasta

#klf
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.1.2_mammalia_dbd_fasta.fasta

#sp1
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/2.3.1.1_mammalia_dbd_fasta.fasta

#cebp
wget http://tfclass.bioinf.med.uni-goettingen.de/suppl/1.1.8.1_mammalia_dbd_fasta.fasta
echo '>Mus_musculus_CEBPD-annot_ma-DBD_ma' >> 1.1.8.1_mammalia_dbd_fasta.fasta
echo 'RQRRERNNIAVRKSRDKAKRRNQEMQQKLVELSAENEKLHQRVEQLTRDLAGLRQFFK' >> 1.1.8.1_mammalia_dbd_fasta.fasta

#concatenate
cat 1.1.1*mammalia_dbd_fasta.fasta 1.1.2*mammalia_dbd_fasta.fasta 1.1.3*mammalia_dbd_fasta.fasta 1.1.4*mammalia_dbd_fasta.fasta | grep -A 1 ">Mus_musculus_" | grep -v '\--' > AP1_dbd_class.fasta
cat 2.1.1.1_mammalia_dbd_fasta.fasta | grep -A 1 ">Mus_musculus_" | grep -v '\--' > GR_dbd_class.fasta
cat 2.3.1.2_mammalia_dbd_fasta.fasta | grep -A 1 ">Mus_musculus_" | grep -v '\--' > KLF_dbd_class.fasta
cat 2.3.3.16_mammalia_dbd_fasta.fasta 1.2.3.2_mammalia_dbd_fasta.fasta  | grep -A 1 ">Mus_musculus_" | grep -v '\--' > TWIST-ZBTB_dbd_class.fasta
cat 2.3.1.1_mammalia_dbd_fasta.fasta | grep -A 1 ">Mus_musculus_" | grep -v '\--' > SP_dbd_class.fasta
cat 1.1.8.1_mammalia_dbd_fasta.fasta | grep -A 1 ">Mus_musculus_" | grep -v '\--' > CEBP_dbd_class.fasta

cat AP1_dbd_class.fasta GR_dbd_class.fasta KLF_dbd_class.fasta SP_dbd_class.fasta TWIST-ZBTB_dbd_class.fasta CEBP_dbd_class.fasta > all_dbd_classes.fasta

/Users/guertinlab/Downloads/fasta-36.3.8h/bin/ssearch36 -s MD40 -m 8CBl all_dbd_classes.fasta all_dbd_classes.fasta > test_output.txt

grep '>' all_dbd_classes.fasta > order_file.txt


#pswms
ceqlogo -i KLF_meme.txt -m 1 -o KLF_meme.eps
ceqlogo -i KLF_meme.txt -m 1 -o KLF_meme_rc.eps -r

ceqlogo -i SP_meme.txt -m 1 -o SP_meme.eps
ceqlogo -i SP_meme.txt -m 1 -o SP_meme_rc.eps -r

ceqlogo -i PSWM_family_1_meme.txt -m 1 -o AP1_meme.eps
ceqlogo -i PSWM_family_1_meme.txt -m 1 -o AP1_meme_rc.eps -r

ceqlogo -i PSWM_family_2_meme.txt -m 1 -o bHLH_meme.eps
ceqlogo -i PSWM_family_2_meme.txt -m 1 -o bHLH_meme_rc.eps -r

ceqlogo -i PSWM_family_3_meme.txt -m 1 -o GR_meme.eps
ceqlogo -i PSWM_family_3_meme.txt -m 1 -o GR_meme_rc.eps -r

ceqlogo -i PSWM_family_6_meme.txt -m 1 -o ETS_meme.eps
ceqlogo -i PSWM_family_6_meme.txt -m 1 -o ETS_meme.eps

ceqlogo -i PSWM_family_7_meme.txt -m 1 -o ZBTB33-GFX_meme.eps
ceqlogo -i PSWM_family_7_meme.txt -m 1 -o ZBTB33-GFX_meme.eps

ceqlogo -i PSWM_family_8_meme.txt -m 1 -o Maz_meme.eps
ceqlogo -i PSWM_family_8_meme.txt -m 1 -o Maz_meme.eps

ceqlogo -i PSWM_family_9_meme.txt -m 1 -o NFY_meme.eps
ceqlogo -i PSWM_family_9_meme.txt -m 1 -o NFY_meme_rc.eps -r

ceqlogo -i PSWM_family_10_meme.txt -m 1 -o NRF_meme.eps
ceqlogo -i PSWM_family_10_meme.txt -m 1 -o NRF_meme_rc.eps -r

ceqlogo -i PSWM_family_11_meme.txt -m 1 -o Olig2_meme.eps
ceqlogo -i PSWM_family_11_meme.txt -m 1 -o Olig2_meme_rc.eps -r

ceqlogo -i PSWM_family_12_meme.txt -m 1 -o STAT_meme.eps
ceqlogo -i PSWM_family_12_meme.txt -m 1 -o STAT_meme_rc.eps -r

ceqlogo -i PSWM_family_13_meme.txt -m 1 -o TFAP2A_meme.eps
ceqlogo -i PSWM_family_13_meme.txt -m 1 -o TFAP2A_meme_rc.eps -r

ceqlogo -i PSWM_family_14_meme.txt -m 1 -o TWIST_meme.eps
ceqlogo -i PSWM_family_14_meme.txt -m 1 -o TWIST_meme_rc.eps -r
