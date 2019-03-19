#!/usr/bin/python3
import os,sys,time, urllib.request, tarfile, json
import xml.etree.ElementTree as ET, logging
from helper import *

print("\nAll STDOUT & STDERR are logged in log.txt file\n")
stdout_logger = logging.getLogger('STDOUT')
sys.stdout = StreamToLogger(stdout_logger, logging.INFO)

stderr_logger = logging.getLogger('STDERR')
sys.stderr = StreamToLogger(stderr_logger, logging.ERROR)

if len(sys.argv) != 2:
    print("python3 <script> <jsonfile(file with query terms)>\n Exiting!!!!\n")
    exit(0)

if not os.path.exists(sys.argv[1]):
    print("Input file not found\nExiting!!\n")
    exit(0)

def parse_familyxml(file,gse,gsmid_processed):#,gsmfile):
    outfile = open("sc_protocol_rnaseq.tsv",'a')
    outfile1 = open("sc_protocol_complete.tsv","a")
    root  = ET.parse(file).getroot()
    gsmid_new = []
    uniq_protocol = []
    minxml_prefix = "{http://www.ncbi.nlm.nih.gov/geo/info/MINiML}"
    for child in root:
        #print(child.tag)
        if child.tag == minxml_prefix+"Sample":
            gsm = child.attrib["iid"].strip()
            if gsm in gsmid_processed:
                continue
            else:
                gsmid_processed.append(gsm)
                gsmid_new.append(gsm)
            charecteristic_tag= []
            instru=[]
            title = gsmid = source = organism = growth = molecule = \
                    protocol = dataprocess = libstrat = libsource = \
                    libsel =  biosample = srx = char_tag = instru_str =""
            for grandchild in child:
                #print("\tgrandchild:",grandchild.tag)
                #print("char:",char_tag)
                #print("isntru:",instru_str)
                if grandchild.tag == minxml_prefix+"Title":
                    title = grandchild.text.strip()
                elif grandchild.tag == minxml_prefix+"Accession":
                    gsmid = grandchild.text.strip()
                elif grandchild.tag == minxml_prefix+"Channel":
                    for i in grandchild:
                        #print("\t\ti:",i.tag)
                        if i.tag == minxml_prefix+"Source":
                            source = i.text.strip()
                        elif i.tag == minxml_prefix+"Organism":
                            organism = i.text.strip()
                        elif i.tag == minxml_prefix+"Characteristics":
                            try:
                                charecteristic_tag.append(i.attrib["tag"].strip()+":"+i.text.strip())
                            except KeyError:
                                 charecteristic_tag.append("NA"+":"+i.text.strip())
                        elif i.tag == minxml_prefix+"Growth-Protocol":
                            growth = i.text.strip().replace('\n', '').replace('\r','')
                        elif i.tag == minxml_prefix+"Molecule":
                            molecule = i.text.strip()
                        elif i.tag == minxml_prefix+"Extract-Protocol":
                            protocol = i.text.strip().replace('\n', '').replace('\r','')
                            if protocol not in uniq_protocol:
                                uniq_protocol.append(protocol)
                elif grandchild.tag == minxml_prefix+"Data-Processing":
                    dataprocess = grandchild.text.strip().replace('\n', '').replace('\r','')
                elif grandchild.tag == minxml_prefix+"Library-Strategy":
                    libstrat = grandchild.text.strip()
                elif grandchild.tag == minxml_prefix+ "Library-Source":
                    libsource = grandchild.text.strip()
                elif grandchild.tag ==  minxml_prefix+"Library-Selection":
                    libsel = grandchild.text.strip()
                elif grandchild.tag == minxml_prefix+  "Instrument-Model":
                    for i in grandchild:
                        instru.append(i.text.strip())
                elif grandchild.tag == minxml_prefix+ "Relation" and \
                        grandchild.attrib["type"] == "Biosample":
                    biosample = grandchild.attrib["target"].strip()
                elif grandchild.tag == minxml_prefix+"Relation" and \
                        grandchild.attrib["type"] == "SRA":
                    srx  = grandchild.attrib["target"].strip()
            char_tag = ";".join(charecteristic_tag)
            instru_str = ";".join(instru)
            out_line =gse+"\t"+gsmid+"\t"+title+"\t"+source+"\t"+organism+"\t"+ \
                              libstrat +"\t"+libsource +"\t"+libsel+"\t"+ \
                              char_tag+"\t"+instru_str+"\t"+ \
                             growth+"\t"+\
                             molecule+"\t"+\
                            protocol+"\t"+\
                            dataprocess+"\t"+\
                              biosample+"\t"+srx+"\t"+"\n"
            if libstrat == "RNA-Seq":
                outfile.write(out_line)
            outfile1.write(out_line)
    outfile.close()
    outfile1.close()
    print("\tparsing family tree done",file=sys.stdout)
    return ";".join(uniq_protocol), gsmid_new


def unzip_tar(tgzfile,path):
    if tarfile.is_tarfile(tgzfile):
        if not os.path.exists(tgzfile.strip(".tgz")):
            tar = tarfile.open(tgzfile,"r:gz")
            tar.extractall(path)
            tar.close()
        print("\tparsing "+tgzfile.strip(".tgz"),file=sys.stdout)
        return(True)
    else:
        print("\t "+ tgzfile + "not a tar file",file=sys.stdout)
        return(False)


def parse_protocol(ftplink):
    out = open("noGSE_protocolxml.txt",'a')
    gsmid_processed = []
    #gsmfile = "completed_gsm.txt"
    exist = False
    if os.path.exists(gsmfile):
        print("\tFound processed gsmid file",file=sys.stdout)
        gsmin = open(gsmfile)
        for line in gsmin:
            if line!="" and line != "\n":
                gsmid_processed.append(line.strip())
        gsmin.close()
        exist = True


    basename = ftplink.split("/")[-2]
    protocol_link = ftplink + "miniml/" +basename+"_family.xml.tgz"
    saveDir = "protocol/"
    file = saveDir+os.path.basename(protocol_link)
    #print(file)
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    if not os.path.exists(file):
        try:
            urllib.request.urlretrieve(protocol_link,file)
            print("\tDownloaded protocol file",file=sys.stdout)
        except (ValueError, urllib.error.URLError ):
            print("=====================",file=sys.stdout)
            print("\t"+basename+" file not found",file=sys.stdout)
            print(ftplink)
            print("=====================\n",file=sys.stdout)
            out.write(basename+"\t"+ftplink+"\n")
            out.close()
            return("",[])
    status = unzip_tar(file,saveDir)
    if(status or exist == True):
        print("\tParsing family xml",file=sys.stdout)
        proto, l = parse_familyxml(file.strip(".tgz"),basename,gsmid_processed)#,gsmfile)
        return(proto,l)
    else:
        print("\t"+basename+" file not found",file=sys.stdout)
        out.write(basename+"\t"+ftplink+"\n")
        out.close()
        return("",[])
    out.close()


def find_library_prep(title,abstract,protocol):
    library = {
        "smart-seq":"smart-seq",
        "smartseq": "smart-seq",
        "smartseq2": "smart-seq2",
        "smart-seq2": "smart-seq2",
        "10x" : "10x",
        "chromium" : "10x",
        "drop-seq": "drop-seq",
        "indrop" : "indrop",
        "droplet": "droplet",
        "microwell": "microwell",
        "cel-seq": "cel-seq",
        "cel-seq2": "cel-seq2",
        "quartz-seq": "quartz-seq",
        "mars-seq" : "mars-seq",
        "scrb-seq" : "scrb-seq",
        "strt" : "strt-seq",
        "fluidigm": "fluidigm"

    }
    library_prep = []
    #equipment = []
    for i in library.keys():
        if title.lower().find(i)!=-1 or abstract.lower().find(i)!=-1 or \
           protocol.lower().find(i) != -1:
            library_prep.append(library[i])
    if len(library_prep)==0:
        library_prep.append("NA")
    return(";".join(library_prep))


def parse_esummary(root):
    outputfile="Geo_sc_Datasets.tsv"
    outputfile1="Geo_10x_dataset.tsv"
    if os.path.exists(uidfile):
        outfile = open(outputfile,"a")
    else:
        outfile = open(outputfile,"w")
        outfile.write("uid\tgse\tsrp\tspecies\tlib_prep\tnoofsamples\ttitle"+
                      "\tftplink\t"+
                      "targetftplink\tpubmedid"+"\thumandata\tmousedata\t"+
                      "otherdata\tgplid\tabstract\n")
    if os.path.exists(uidfile):
        outfile1 = open(outputfile1,"a")
    else:
        outfile1 = open(outputfile1,"w")
        outfile1.write("uid\tgse\tsrp\tspecies\tlib_prep\tnoofsamples\ttitle"+
                       "\tftplink\t"+
                      "targetftplink\tpubmedid"+"\thumandata\tmousedata\t"+
                      "otherdata\tgplid\tabstract\n")
    outfile2 = open(uidfile,'a')

    for child in root:
        gsmid_new =[]
        uid = gse= title = abstract = species = gplid = srp = ftplink = \
                targetftplink = noofsamples = pubmedid = uniq_protocol = ""
        human = "False"
        mouse = "False"
        other = "False"

        for grandchild in child:
            if grandchild.tag == "Id":
                uid = grandchild.text.strip()
                print("\n===========\n"+ time.strftime("%d-%m-%Y %H:%M:%S",
                                         time.localtime()) + "\n",file=sys.stdout)
                print(uid,file=sys.stdout)
                continue
            elif grandchild.tag == "Item" and grandchild.attrib["Name"] == "Accession":
                gse=grandchild.text
                print(gse,file=sys.stdout)
                continue
            elif grandchild.tag == "Item" and grandchild.attrib["Name"] == "title":
                title = grandchild.text.strip().replace('\n','').replace('\r','')
                continue
            elif grandchild.tag == "Item" and grandchild.attrib["Name"] == "summary":
                abstract = grandchild.text.strip().replace('\n','').replace('\r','')
                continue
            elif grandchild.tag == "Item" and grandchild.attrib["Name"] == "taxon":
                #print("taxon")
                species=grandchild.text.strip()
                if species.lower().find("homo sapiens") != -1:
                    human = "Yes"
                if species.lower().find("mus musculus") != -1:
                    mouse = "Yes"
                if species.lower().find("homo sapiens") == -1 and  \
                   species.lower().find("mus musculus") == -1:
                    other = "Yes"
                continue
            elif grandchild.tag == "Item" and grandchild.attrib["Name"] == "GPL":
                gplid=grandchild.text
                continue
            elif grandchild.tag == "Item" and grandchild.attrib["Name"] == "ExtRelations":
                for i in grandchild:
                    for j in i:
                        if j.tag == "Item" and j.attrib["Name"] == "TargetObject":
                            srp= j.text
                        elif j.tag == "Item" and j.attrib["Name"] == "TargetFTPLink":
                            targetftplink = j.text.strip()
                continue
            elif grandchild.tag == "Item" and grandchild.attrib["Name"] == "n_samples":
                noofsamples=grandchild.text
                continue
            elif grandchild.tag == "Item" and grandchild.attrib["Name"] == "PubMedIds":
                pubmedlist=[]
                for i in grandchild:
                    if i.tag == "Item" and i.attrib["Name"] =="int":
                        pubmedlist.append(i.text)
                if len(pubmedlist) > 1:
                    for i in pubmedlist:
                        pubmedid = pubmedid + i + ";"
                elif len(pubmedlist) == 1:
                    pubmedid = pubmedlist[0]
                else:
                    pubmedid = ""
                continue
            elif grandchild.tag == "Item" and grandchild.attrib["Name"] == "FTPLink":
                ftplink=grandchild.text.strip()
                #print(ftplink)
                # Download and parse protocol files
                uniq_protocol, gsmid_new = parse_protocol(ftplink)
                #print(uniq_protocol)
                #print(len(gsmid_new))

                #print(pubmedid)

        lib_prep = find_library_prep(title,abstract,uniq_protocol)
        print("\tFinding which lib prep done",file=sys.stdout)
        if lib_prep.find("10x") != -1:
            outfile1.write(str(uid)+"\t"+str(gse)+"\t"+str(srp)+"\t"+species+"\t"+
                      lib_prep+"\t"+str(noofsamples)+"\t"+title+
                      "\t"+str(ftplink)+"\t"+
                      targetftplink+"\t"+str(pubmedid)+"\t"+str(human)+"\t"+str(mouse)+"\t"+
                      str(other)+"\t"+str(gplid)+
                           "\t"+abstract+"\n")
        #print(uid,gse,species,gplid)
        outfile.write(str(uid)+"\t"+str(gse)+"\t"+str(srp)+"\t"+species+"\t"+
                      lib_prep+"\t"+str(noofsamples)+"\t"+title+
                      "\t"+str(ftplink)+"\t"+
                      targetftplink+"\t"+str(pubmedid)+"\t"+str(human)+"\t"+str(mouse)+"\t"+
                      str(other)+"\t"+str(gplid)+
                      "\t"+abstract+"\n")
        outfile2.write(uid+"\n")
        if len(gsmid_new) > 0:
            gsmout = open(gsmfile,'a')
            for i in gsmid_new:
                gsmout.write(i+"\n")
            gsmout.close()
        else:
            outhandle = open("noGSE_forGSMS.txt",'a')
            outhandle.write(gse+"\n")
            outhandle.close()
            print(gse+" No new GSMs",file=sys.stdout)
    print("\n===========",file=sys.stdout)

    outfile.close()
    outfile2.close()
    outfile1.close()

def esearch_ncbi(esearch):
# function to search gds based on query term & return xml

    global uid_list
    global prefix
    try:
        esearch = urllib.request.urlopen(prefix+esearch)
    except urllib.error.URLError:
        #ResponseData = e.read().decode("utf-8", 'ignore')
        print("URL Error occured! Check connection\n",file=sys.stderr)
        return
    except urllib.error.HTTPerror():
        print("HTTP Error occured! Check URL\n",file=sys.stderr)
        return
    esearch.xml = esearch.read().decode('utf-8')
    root  = ET.fromstring(esearch.xml)

#    for child in root:
#        if child.tag=="Count": count = int(child.text.strip())
#        if child.tag=="RetMax": retmax = int(child.text.strip())
#        if child.tag=="RetStart": retstart = int(child.text.strip())
#        if child.tag=="QueryKey": query_key = child.text.strip()
#        if child.tag=="WebEnv": webenv =child.text.strip()
    for id in root.iter("Id"):
        if id.text.strip() not in uid_list:
            uid_list.append(id.text.strip())


################# main funtion###############
if __name__ == "__main__":

    #print ("Start : %s" % time.ctime(),file=sys.stdout)
    preface = """###############################################
START: {time}
"""
    print(preface.format(time = time.ctime()),file=sys.stdout)


    uid_list = []
    uid_processed = []
    uid_toprocess = []

    uidfile="uid_finished.txt"
    gsmfile="completed_gsm.txt"

    if os.path.exists(uidfile):
        print("Existing uid file found",file=sys.stdout)
        uid_finished = open(uidfile)
        for line in uid_finished:
            if line!="" and line != "\n":
                uid_processed.append(line.strip())
    else:
        print("No processed uids",file=sys.stdout)
# Note: if you are searching for chip-seq/atac-seq some other non expression
# the esearch suffix has to be changed

    prefix='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    esearch_gds='esearch.fcgi?db=gds&term='
    esearch_suffix ='+AND+expression+profiling+by+high+throughput+sequencing[DataSet%20Type]&retmax=5000'

############################# Query terms #########################
###################################################################
## queries to search

    #query_terms = [
    #'"single+cell+RNA%2dseq"',
    #    '"single+cell+RNAseq"','"scRNAseq"','"single-cell+RNAseq"',
    #    '"singlecellRNAseq"','singlecell+RNAseq',"single_cell_rnaseq",
    #    "singlecell_rnaseq"

    #]
    #with open(sys.argv[1],"w") as filehandle:
    #    json.dump(query_terms,filehandle)
    #exit(0)
    with open(sys.argv[1],"r") as filehandle:
        query_terms = json.load(filehandle)
    
    i=1
    for term in query_terms:
        print(str(i)+":searching  query:"+term,file=sys.stdout,sep="")
        esearch_ncbi(esearch_gds+term+esearch_suffix)
        print("After "+str(i)+" search hits  " + str(len(uid_list)),file=sys.stdout,sep="")
        time.sleep(5)
        i=i+1
    #time.sleep(5)

#################################################################

    uid_toprocess = list(set(uid_list) - set(uid_processed))
    #uid_toprocess = uid_toprocess[0:35]

#print("Count:",str(count))
#print("Retmax:",str(retmax))
#print("QueryKey:",query_key)
#print("WebEnv:",webenv)
    print("No of uids in search:"+str(len(uid_list)))
    print("No of uids already processed:" + str(len(uid_processed)))
    print("No of uids to process:" + str(len(uid_toprocess)))

    if len(uid_toprocess)==0:
        print("Nothing new to process",file=sys.stdout)
        exit()

# if for some reason a uid fails -> rerun -> it will resume:
# will download only uids that are not completed
# if the same uid fails again, you might want to add to add it to
# finished uid list

    start=0
    step=30
    print("To process:",str(len(uid_toprocess)))

    while(start < len(uid_toprocess)):
        end = start +step
        if end > len(uid_toprocess):
            end= len(uid_toprocess)
        print("start="+str(start)+" end="+str(end),file=sys.stdout)
        queryids = ','.join(uid_toprocess[start:end])
        url = prefix + "esummary.fcgi?db=gds&id=" + queryids
        esummary =  urllib.request.urlopen(url)
        root  =   ET.fromstring(esummary.read().decode('utf-8'))
        parse_esummary(root)
        start = start + step
        time.sleep(5)
        #exit(0)


    print ("End : %s" % time.ctime(),file = sys.stdout)
