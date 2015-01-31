import sys

ped_in = open(sys.argv[1])
dat = open(sys.argv[2],"r")
number = int(sys.argv[3])

#store family information and maternal and paternal genotypes separately
family=[]
all_maternal=[]
all_paternal=[]
for i in ped_in:
    s=i.rstrip("\n").split()
    count=0
    maternal=[]
    paternal=[]
    if s[0]!="end":
        family.append(s[0:5])
    gen_count=0
    loop_count=0
    for j in s:
        count=count+1
        if count>5:
            gen_count=gen_count+1
            if gen_count==loop_count*(number+1)+1:
                loop_count=loop_count+1
            else:
                gen=j.split("/")
                if len(gen)>1:
                    maternal.append(gen[0])
                    paternal.append(gen[1])
    all_maternal.append(maternal)
    all_paternal.append(paternal)

#concatenate maternal alleles
this_grow=[]
this=""
maternal_grow=[]
for i in all_maternal:
    count=0
    for j in i:
        count=count+1
        if count==1:
            this=str(j)
        if count>1 and count<(number+1):
            this=this+"."+str(j)
        if count==(number+1):
            this_grow.append(this)
            count=1
            this=str(j)
    if len(this_grow)!=0:
        maternal_grow.append(this_grow)
    this_grow=[]


#concatenate paternal alleles
this_grow=[]
this=""
paternal_grow=[]
for i in all_paternal:
    count=0
    for j in i:
        count=count+1
        if count==1:
            this=str(j)
        if count>1 and count<(number+1):
            this=this+"."+str(j)
        if count==(number+1):
            this_grow.append(this)
            count=1
            this=str(j)
    if len(this_grow)!=0:
        paternal_grow.append(this_grow)
    this_grow=[]

#store founders
count=-1
maternal_founders=[]
paternal_founders=[]
founders=[]
id=[]
for i in family:
    id.append(i[1])
    count=count+1
    if i[2]=="0" and i[3]=="0":
        paternal_founders.append(paternal_grow[count])
        maternal_founders.append(maternal_grow[count])
        founders.append(i[1])


#store SNP names

dat_store=dict([])
count=0
for i in dat:
    s=i.strip("\n").split()
    if s[1][0] != "M":
        count=count+1
        dat_store[count]=s[1]

allele_store=dict([])
#store maternal genotypes
count=0
for i in zip(*maternal_founders):
    count=count+1
    id_count=0
    occurences=[]
    for j in list(i):
        occurences.append(list(i).count(j))
    if max(occurences)==1:
        for j in list(i):
            allele_store[dat_store[count]+"_"+j]="M_"+founders[id_count]
            id_count=id_count+1
    else:
        print "non-unique founder sequence"

#store paternal genotypes
count=0
for i in zip(*paternal_founders):
    count=count+1
    id_count=0
    occurences=[]
    for j in list(i):
        occurences.append(list(i).count(j))
    if max(occurences)==1:
        for j in list(i):
            allele_store[dat_store[count]+"_"+j]="P_"+founders[id_count]
            id_count=id_count+1
    else:
        print "non-unique founder sequence"
        

#output allele dictionary allocations
#output=open("allele_dictionary_out.txt","w")
#for i in allele_store.keys():
#    output.write(i+"\t"+allele_store[i]+"\n")
#output.close()

out_maternal=open(sys.argv[4]+"_maternal.txt","w")
id_out=["SNP"]
count=0
line_out=[]
for i in zip(*maternal_grow):
    count=count+1
    id_count=0
    if count==1:
        for j in list(i):
            id_out.append(id[id_count])
            id_count=id_count+1
        out_maternal.write("\t".join(id_out)+"\n")
    snp=dat_store[count]
    for j in list(i):
        retrieved=allele_store.get(dat_store[count]+"_"+j,"?")
        #print dat_store[count]+"_"+j, retrieved
        line_out.append(retrieved)
    line_out.insert(0,snp)
    out_maternal.write("\t".join(line_out)+"\n")
    line_out=[]
out_maternal.close()
        
out_paternal=open(sys.argv[4]+"_paternal.txt","w")
id_out=["SNP"]
count=0
line_out=[]
for i in zip(*paternal_grow):
    count=count+1
    id_count=0
    if count==1:
        for j in list(i):
            id_out.append(id[id_count])
            id_count=id_count+1
        out_paternal.write("\t".join(id_out)+"\n")
    snp=dat_store[count]
    for j in list(i):
        retrieved=allele_store.get(dat_store[count]+"_"+j,"?")
        line_out.append(retrieved)
    line_out.insert(0,snp)
    out_paternal.write("\t".join(line_out)+"\n")
    line_out=[]
out_paternal.close()
        
