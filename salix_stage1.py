import sys, time

print
print "\t##################################"
print "\t####  Haplotype Flow Drop v1  ####"
print "\t####       A. P. Levine       ####"
print "\t##################################"
print

print "\t" + time.strftime("%A %d/%m/%Y %H:%M:%S", time.localtime())

profile = open(sys.argv[1],"r")
cm_lengths = open(sys.argv[2],"r")

spacing = float(sys.argv[3])

print

print "\tPedigree: "+ sys.argv[1]

print

print "Stage 1 - cM files:"
print "\t" + "Preparing cM MAP with " + str(spacing) + " cM spacing ..."

mapout = open("genome.map","w")
map_store = []
count = 0
for i in cm_lengths:
	s=i.rstrip("\n").split()
	chromosome = s[0]
	cm = round(float(s[1]))+1
	now=0-spacing
	while now<cm:
		count=count+1
		now = now+spacing
		marker = "SNP_"+str(count)
		position = str(now)
		mapout.write(str(chromosome)+"\t"+marker+"\t"+position+"\n")
		map_store.append([chromosome,marker,position])
mapout.close()
map_count = len(map_store)	

print "\tPreparing cM PED file ..."

#Prepare cM PED file
pedfile = open("genome.ped","w")
for i in profile:
	s=i.rstrip("\n").split()
	pedfile.write(" ".join(s) + ' 1 1'*map_count+"\n")
pedfile.close()

print "\tPreparing cM FRQ file ..."

#Prepare cM FREQ file
freqfile = open("genome.frq","w")
for i in map_store:
	freqfile.write("M" + "\t" + i[1] +"\n")
	freqfile.write("A" + "\t" + "1" + "\t" + "0.5" +"\n")
	freqfile.write("A" + "\t" + "2" + "\t" + "0.5" +"\n")
freqfile.close()

print "\tPreparing cM DAT file ..."

#Prepare cM DAT file
datfile = open("genome.dat","w")
datfile.write("A\tAnalysis\n")
for i in map_store:
	datfile.write("M" + "\t" + i[1] +"\n")
datfile.close()

cM_next=0
number_alleles=int(sys.argv[4])
out_number=int(sys.argv[5])

print
print "Stage 2 - marker tagging:"
print "\tIncoporating " + str(out_number) + " tag markers each with "+ str(number_alleles) + " alleles at "+ str(cM_next) +" cM spacing for each SNP ..."

#Parameters
pedfile=open("genome.ped","r")
datfile=open("genome.dat","r")
mapfile=open("genome.map","r")
frqfile=open("genome.frq","r")

out_name="for_simulation"

#Map file
print "\tPreparing " + out_name +".map ..."
new_mapfile = open(out_name+".map","w")
new_tags = open(out_name+".tag","w")
count=0
for i in mapfile:
    s=i.rstrip("\n").split()
    previous=s
    new_mapfile.write(" ".join(s)+"\n")
    loop_count=0
    grow=[]
    grow.append(s[1])
    for j in range(out_number):
        loop_count=loop_count+1
        count=count+1
        new_mapfile.write(previous[0]+" M"+str(count)+" "+str(float(previous[2])+cM_next*loop_count)+"\n")
        grow.append("M"+str(count))
    new_tags.write(" ".join(grow)+"\n")
new_mapfile.close()
new_tags.close()

#Frequency file
print "\tPreparing " + out_name +".frq ..."
new_freqfile = open(out_name+".frq","w")
each_freq=str(1.0/number_alleles)
first=0
top=1
count=0
go=0
for i in frqfile:
    s=i.rstrip("\n").split()
    if first==0:
        new_freqfile.write(" ".join(s)+"\n")
    if first==1:
        if s[0]!="M" and top==1:
            new_freqfile.write(" ".join(s)+"\n")
        if s[0]=="M":
            go=0
            top=0
            for j in range(out_number):
                count=count+1
                new_freqfile.write("M M"+str(count)+"\n")
                for k in range(0,number_alleles):
                    new_freqfile.write("A "+str(k+1)+" "+each_freq+"\n")
            new_freqfile.write(" ".join(s)+"\n")
            go=1
        if go==1 and s[0]!="M":
            new_freqfile.write(" ".join(s)+"\n")
    first=1
for j in range(out_number):
    count=count+1
    new_freqfile.write("M M"+str(count)+"\n")
    for k in range(0,number_alleles):
        new_freqfile.write("A "+str(k+1)+" "+each_freq+"\n")
new_freqfile.close()

#Dat file
print "\tPreparing " + out_name +".dat ..."
new_datfile = open(out_name+".dat","w")
count=0
first=0
for i in datfile:
    s=i.rstrip("\n").split()
    if first==0:
        new_datfile.write(" ".join(s)+"\n")
    if first==1:
        new_datfile.write(" ".join(s)+"\n")
        for j in range(out_number):
            count=count+1
            new_datfile.write("M M"+str(count)+"\n")
    first=1
new_datfile.close()

#Ped file
print "\tPreparing " + out_name +".ped ..."
new_pedfile = open(out_name+".ped","w")
for i in pedfile:
    s=i.rstrip("\n").split()
    count=0
    line=[]
    for j in s:
        count=count+1
        if count <7:
            line.append(j)
        if count >6:
            line.append(j)
            if count%2==0:
                for j in range(out_number):
                    line.append('1')
                    line.append('1')
    new_pedfile.write(" ".join(line)+"\n")
new_pedfile.close()

print
print "Finished at " + time.strftime("%A %d/%m/%Y %H:%M:%S", time.localtime())
