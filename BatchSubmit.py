import os

#########################
#types to submit
#########################
types={}
types["dijet"] = ["0",5,100]
types["wz"]    = ["1",5,100]
types["ttbar"] = ["2",5,100]

######################
#Clean out batch directory
######################
os.system("rm -r BatchOutput")
os.system("mkdir BatchOutput")

transferdir = os.getcwd()+"/BatchOutput"

######################
#Submit all jobs
######################
for type in sorted(types.keys()):

    typeflag     = types[type][0]
    njobs        = types[type][1]
    eventsperjob = str(types[type][2])

    for i in range(njobs):
        print "\n********************"
        print "Submitting job: ",i
    
        #change to bsub on TeVCluster?
        command = "source BatchSubmit.sh "+typeflag+" "+type+"_"+str(i)+" "+eventsperjob+" "+transferdir
        print command
        os.system(command)