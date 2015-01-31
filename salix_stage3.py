import sys
import random

#hapdrop_v5_unphased

#1	pedigree.pro
#2	positions_16_mapped.txt
#3	sim_flow_1_maternal.txt
#4	sim_flow_1_paternal.txt
#5	genotypes16.txt
#6	output

print "Storing pedigree information..."
print "Identifying founders..."
input_pedigree=open(sys.argv[1],"r")
pedigree=[]
pedigree_id=[]
founders = []
affected = []
for i in input_pedigree:
    s=i.rstrip("\n").split()
    pedigree.append(s)
    pedigree_id.append(s[1])
    if s[2] == "0" and s[3] =="0": #founders
	founders.append(s[1])
    if s[5] == "2":
	affected.append(s[1])

#founders = []
#affected = []
#store=dict([])
#for i in open(sys.argv[1],"r"): #load pedigree (e.g. pedigree.pro)
#	s=i.rstrip("\n").split()
#	store[s[1]]=[s[2],s[3]]
#	if s[2] == "0" and s[3] =="0": #founders
#		founders.append(s[1])
#	if s[5] == "2": #affected
#		affected.append(s[1])
#affected_ancestors = []
#for i in affected:
#	i_founders = store.get(i,"NA")
#	for j in i_founders:
#		if j!="0":
#			retrieved = store.get(j,"NA")
#			if j!="NA":
#				i_founders.append(retrieved[0])
#				i_founders.append(retrieved[1])
#	for j in i_founders:
#		if j not in affected_ancestors and j != "0":
#			affected_ancestors.append(j)
#aff_founders = []
#for i in founders:
#	if i in affected_ancestors:
#		aff_founders.append(i)




print "Preparing SNP mapping dictionary..."
input_found=open(sys.argv[2],"r") #load file mapping each variant to SNP (e.g. positions_16_mapped.txt)
found=dict([])
for i in input_found:
    s=i.rstrip("\n").split()
    chr=s[0]
    pos=s[1]
    snp=s[3]
    found[chr+"."+pos]=snp #convert from chr.pos to SNP ID

print "Storing haplotype flow file 1 data..."
input_sim1=open(sys.argv[3],"r") #haplotype file 1 (e.g. sim_flow_1_maternal.txt), it does not matter whether maternal or paternal is loaded first
first=0
header=[] #individuals in simulated haplotype file (ordered)
sim1=dict([])
for i in input_sim1:
    s=i.rstrip("\n").split()
    if first==0:
        for j in s:
                header.append(j)
        header.remove("SNP")
    if first==1:
        sim1[s[0]]=s[1:]
    first=1

#For each SNP, store the founder source of each persons allele in the order in which it appears

print "Storing haplotype flow file 2 data..."
input_sim2=open(sys.argv[4],"r") #haplotype file 2 (e.g. sim_flow_1_paternal.txt)
first=0
sim2=dict([])
for i in input_sim2:
    s=i.rstrip("\n").split()
    if first==1:
        sim2[s[0]]=s[1:]
    first=1

#print "Storing sequenced individuals..."
#sequenced=[]
#input_sequenced=open(sys.argv[5],"r")
#for i in input_sequenced:
#    s=i.rstrip("\n").split()
#    sequenced.append(s[0])

def divide(numerator, denominator):
    if denominator==0:
        return 0
    else:
        return float(numerator)/denominator
   

print "Randomising founder genotypes, identifying, mapping..."
input_genotypes=open(sys.argv[5],"r") #genotype data of controls (e.g. genotypes16.txt)
first=0
grow_ped=[]
grow_txt=[]
grow_txt_sequenced=[]
keep=[]
keep_sequenced=[]
for i in input_genotypes:
    s=i.rstrip("\n").split("\t")
    allocated_d=dict([])
    if first==0: #header row
        #number=len(s)-5 #number of control individuals
        #controls=range(5,number+5) #column numbers for controls
	controls=range(5,len(s))
        allocated=random.sample(controls,len(founders)) #random subset of controls
	allocated_id=[]
	for j in allocated:
		allocated_id.append(s[j])
        #aff_founders_pos=[]
        #for j in aff_founders: #for each of the founders of affected individuals...
        #    aff_founders_pos.append(founders.index(j)) #append the position of that individual in the list of all founders
    if first==1: #after header
        chr=s[0]
        pos=s[1]
	identifier=s[2]
	reference=s[3]
	alternate=s[4]
        snp=found.get((chr+"."+pos),"NA")
        if snp=="NA":
            print "SNP not found:", chr, pos
            #Chr/pos not found in SNP mapping dictionary, this should not happen
            pass
        else: #Chr/pos mapped to SNP succesfully
            foundersg=[s[j] for j in allocated] #genotypes of all founders based on random allocation of controls to founders
            #foundersg_aff=[foundersg[k] for k in aff_founders_pos] #genotypes of founders of affected individuals
            count00 = foundersg.count("0/0") + foundersg.count("0|0") #wild type
            countNA = foundersg.count("./.") #missing
            if count00 + countNA == len(founders):
                pass
            else: #only proceed if the there is at least one occurence of the variant in the founders
                maternal=dict([])
                paternal=dict([])
                out_ped=[]
                out_txt=[]
                this_keep=[]
                for j in range(len(founders)): #for each founder
                    name=founders[j] #name of the founder

		    	#to remove...
		    #For genotypes that are unphased...
                    genotype=foundersg[j].split("/")
                    mp=random.sample([0,1],2)
                    m=genotype[mp[0]]
                    maternal["M_"+name]=m
                    p=genotype[mp[1]]
                    paternal["P_"+name]=p

		    #genotype=foundersg[j].split("|") #genotype for this founder
		    #m = genotype[0] #this founder's maternal allele
		    #maternal["M_"+name]=m
                    #p = genotype[1] #this founder's paternal allele
                    #paternal["P_"+name]=p
                count=0
                this_keep.append(chr)
                this_keep.append(pos)
                this_keep.append(snp)
		this_keep.append(identifier)
		this_keep.append(reference)
		this_keep.append(alternate)
                this_first=0
                individual=[]
                bad=0
		#sim1[snp] gives the founder source for each individual for allele 1
		#sim2[snp] ... allele 2
                for j,k in zip(sim1[snp],sim2[snp]): #place the founder source for each allele next to each other
        	            #print j,k
                    this_out_txt="?"
                    this_out_ped="?"
                    this1=maternal.get(j,"NF") #retrieve the actual genotype allele based on the founder allocation
		    #NF: not found
	                    #print "this1 m", this1
                    if this1=="NF": #it could be paternal, not maternal, check for this
                        this1=paternal.get(j,"NF")
        	                #print "this1 p", this1
                    this2=maternal.get(k,"NF")
	     	    #print "this2 m", this2
                    if this2=="NF":
                        this2=paternal.get(k,"NF")
                        #print "this2 p", this2
                    if this1=="NF" or this2=="NF": #founder not found, maternal or paternal, this is a problem!
                        #print sim1[snp]
                        #print sim2[snp]
                        if bad==0:
                            print "not found",chr,pos,snp
                        bad=1
                    if this1=="." or this2=="." or this1=="NA" or this2=="NA": #unknown
                        #this_out_ped="0\t0"
			this_out_ped="0 0" #change to space delimited ped file
                        this_out_txt="NA"
                    if this1=="0" and this2=="0": #wild type
                        #this_out_ped="1\t1"
			this_out_ped="1 1"
                        this_out_txt=0
                    if this1=="1" and this2=="0": #heterozygous
                        #this_out_ped="1\t2"
			this_out_ped="1 2"
                        this_out_txt=1
                    if this1=="0" and this2=="1": #heterozygous
                        #this_out_ped="1\t2"
			this_out_ped="1 2"
                        this_out_txt=1
                    if this1=="1" and this2=="1": #homozygous
                        #this_out_ped="2\t2"
			this_out_ped="2 2"
                        this_out_txt=2
                    if this_first==0:	#if the first time, save the name of the individuals in the order in which they appear
                        individual.append(header[count]) #header from sim file
                        count=count+1
                    #print j,k,this1,this2,this_out_txt
                    out_ped.append(this_out_ped)
                    out_txt.append(this_out_txt)
                if bad==0: #genotype data located succesfully

			#only keep the variant if it is in more than x% of the sequenced individuals
                  #  sequenced_sum=0
                  #  sequenced_g=[]
                  #  for j in sequenced: #for each sequenced individual
                  #      this_sequenced=out_txt[individual.index(j)] #genotype for the sequenced individual
                  #      if this_sequenced!="NA":
                  #          sequenced_sum=sequenced_sum+1 #denominator
                  #          sequenced_g.append(this_sequenced) #genotype for numerator
                  #  frequency=divide(sequenced_sum-sequenced_g.count(0),sequenced_sum) #proportion not wild type
                  #  if frequency>0.5: #if more than 50% of sequenced individuals have the variant...
                  #      grow_txt_sequenced.append(out_txt) #genotypes to keep
                  #      keep_sequenced.append(this_keep) #identity of variant kept


		    #only keep the variant if it is in at least 1 affected
		    to_save = 0
	    	    for j in affected:
			this_affected = out_txt[individual.index(j)]
			if this_affected == 1 or this_affected ==2:
				to_save = 1
		    if to_save == 1:
	                    grow_ped.append(out_ped)
	                    grow_txt.append(out_txt)
	                    keep.append(this_keep)

                    this_first=1

    first=1

pedigree_o=[] #pedigree ordered
for i in individual:
    pedigree_o.append(pedigree[pedigree_id.index(i)])

print "Rotating..."
tgrow_ped=zip(*grow_ped)
tgrow_txt=zip(*grow_txt)
tgrow_txt_sequenced=zip(*grow_txt_sequenced)

print "Saving results..."
#txt file
output=open(sys.argv[6]+".txt","w")
for i,j in zip(pedigree_o,tgrow_txt):
    output.write(" ".join(i)+" "+" ".join(str(k) for k in j)+"\n") #change from tab to space
output.close()
#ped file
output=open(sys.argv[6]+".ped","w")
for i,j in zip(pedigree_o,tgrow_ped):
    output.write(" ".join(i)+" "+" ".join(str(k) for k in j)+"\n") #change from tab to space
output.close()

#only those variants in >x% of sequenced individuals
#output_sequenced=open(sys.argv[6]+"_sequenced","w")
#for i,j in zip(pedigree_o,tgrow_txt_sequenced):
#    output_sequenced.write("\t".join(i)+"\t"+"\t".join(str(k) for k in j)+"\n")
#output_sequenced.close()

#output_snps=open(sys.argv[6]+".variants","w")
output_map=open(sys.argv[6]+".map","w")
for i in keep:
    chr = i[0]
    pos = i[1]
    reference = i[4]
    alternate = i[5]
    if len(reference)>5:
	reference="long"
    if len(alternate)>5:
	alternate="long"
    identifier = i[0]+"_"+i[1]+"_"+i[3]+"_"+reference+"_"+alternate
    #output_snps.write("\t".join(i)+"\n")
    #output_map.write(i[0]+"\t"+i[0]+"."+i[1]+"\t"+i[1]+"\t0\n") #simple map
    output_map.write(chr+"\t"+identifier+"\t"+pos+"\t0\n") #map with identification details
#output_snps.close()
output_map.close()

#output_snps_sequenced=open(sys.argv[6]+"_sequenced_variants","w")
#for i in keep_sequenced:
#    output_snps_sequenced.write("\t".join(i)+"\n")
#output_snps_sequenced.close()

output_allocation=open(sys.argv[6]+".alc","w")
for i in range(len(founders)):
    #output_allocation.write(founders[i]+"\t"+str(allocated[i]-2)+"\n") #what is this -2 for? Perhaps it was for the first 2 columns. Not necessary given that 'controls' from which 'allocated' is defined only includes individuals (not first columns).
    output_allocation.write(founders[i] +"\t"+ str(allocated[i]) +"\t"+allocated_id[i]+"\n")
output_allocation.close()

#print

#print "Files:"
#for i in range(9):
#    print sys.argv[i]



