import sys

def main():
    read_file(sys.argv[1])
    dt = find(sys.argv[1],"dt")
    print("Step size = {} ps".format(dt))
    coupling = find(sys.argv[1], "pcoupltype") + find(sys.argv[1], "pcoupl")
    print(coupling)
    
def find(md_log,param):
    # find value of param in md.log file
    #
    # (eg "dt")
    # in md.log:
    # dt                           = 0.005
    #
    # read md.log file
    file1=open(md_log)
    content=file1.readlines()
    value=0
    param = param + "      "
    for c in content:
        if param in str(c):
            temp=c.split(' = ')[-1]
            try:
                float(temp)
                res=True
            except:
                res=False
            if res:
                value = float(temp)
            else:
                value=temp
    return value
    

def read_file(md_log):
    # read md.log file
    file1=open(md_log)
    content=file1.readlines()

    # using last few lines of any finished md.log file to calculate
    # total simulation time
    cpu_time = float(content[-2].split('  ')[-1])

    # try except statements to deal with the fact if simulations take less
    # than 1 hour (i.e. there is no "h" key character
    hour_per_ns = 0
    try:
        hour_per_ns = float(content[-4].split('h')[0])
    except:
        hour_per_ns =0
    min_per_ns = 0
    try:
        min_per_ns = float(content[-4].split('h')[1].split(':')[0])/60
    except:
        min_per_ns = float(content[-4].split(':')[0])/60
    sec_per_ns = 0
    try:
        sec_per_ns = float(content[-4].split('h')[1].split(':')[1].split('\n')[0])/3600
    except:
        sec_per_ns = float(content[-4].split(':')[-1].split('\n')[0])/3600
        
    print("Total simulation time = {} ns".format(round((hour_per_ns+min_per_ns+sec_per_ns)/cpu_time,3)))


    
main()
