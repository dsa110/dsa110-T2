import numpy as np
import socket 
from T2 import cluster_heimdall, plotting

logger = dsl.DsaSyslogger()
logger.subsystem('software')
logger.app('T2')


def parse_socket(host, port, selectcols=['itime', 'idm', 'ibox', 'ibeam']):
    """ 
    Takes standard MBHeimdall giants socket output and returns full table, classifier inputs and snr tables.
    selectcols will take a subset of the standard MBHeimdall output for cluster. 
    
    host, port: same with heimdall -coincidencer host:port 
    selectcol: list of str.  Select columns for clustering. 
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
    s.bind((host,port))    # assigns the socket with an address
    s.listen(5)             # accept no. of incoming connections
    
    gulp_i = 0 # track heimdall gulp number 
    while True:
        clientsocket, address = s.accept() # stores the socket details in 2 variables
        print(f"Connection from {address} has been established")
        gulp_i += 1 
        print("gulp ", gulp_i)
        
        # read in heimdall socket output  
        ascii_letter = clientsocket.recv(1)           # recieves an alphabet whose ASCII value is the size of the message 
        
        if len(ascii_letter) > 0:
            size = ord(ascii_letter.decode('utf-8'))      # ord() returns the ASCII value of a character
            #candsfile = clientsocket.recv(size)           # recieving the actual msg        
            candsfile = clientsocket.recv(10000000)  # python 3 requires argument here. how to make sure the size is big enough? 
            
            candsfile = candsfile.decode('utf-8')               # decode the bytes msg 
            #print(candsfile)
            
            clientsocket.close()   
            
            # do we want both heimdall.cand and giants.out?  From two sockets?  
            # reading in one gulp at a time. 
            # will modify to continue reading in.  
            tab, data, snrs = cluster_heimdall.parse_candsfile(candsfile)
            
            return tab, data, snrs
        else: 
            print ("this gulp has no heimdall giants output.")



def parse_socket_and_cluster_and_plot(host="127.0.0.1", port=12345, selectcols=['itime', 'idm', 'ibox', 'ibeam'], outputfile="T2_output_socket_", plot=False, plot_dir="socket_"):
    """ 
    continuously read in socket, cluster, dump results to outputfile and save plots.
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM) 
    s.bind((host,port))    # assigns the socket with an address
    s.listen(5)             # accept no. of incoming connections
    
    gulp_i = 0 # track heimdall gulp number 
    
    while True:
        clientsocket, address = s.accept() # stores the socket details in 2 variables
        print(f"Connection from {address} has been established")
        gulp_i += 1 
        print("gulp ", gulp_i)
        
        # read in heimdall socket output  
        candsfile = clientsocket.recv(100000000)  # python 3 requires argument here. how to make sure the size is big enough?         
        
        if len(candsfile) == 0:
            print ("this gulp has no heimdall giants output.")
        
        else: 
            candsfile = candsfile.decode('utf-8')               # decode the bytes msg 
            clientsocket.close()   
           
            print("reading candsfile...")
            print(candsfile.split('\n'))

            tab, data, snrs = cluster_heimdall.parse_candsfile(candsfile)
            
            print("table has", len(tab), "rows")
            
            # T2 cluster 
            clusterer, data_labeled = cluster_heimdall.cluster_data(data, min_cluster_size=10, min_samples=10, metric='euclidean', allow_single_cluster=True, return_clusterer=True)
            clsnr = cluster_heimdall.get_peak(data_labeled, snrs) 
    
            # send T2 cluster results to outputfile
            cluster_heimdall.dump_cluster_results(tab, clsnr, outputfile+str(gulp_i)+".txt", output_cols=['mjds', 'snr', 'ibox', 'dm'])
            
            if plot: 
                plotting.plot_giants(tab, plot_dir=plot_dir+str(gulp_i)+"_") # plot giants      
                plotting.plot_clustered(clusterer, clsnr, snrs, data, tab, cols=['itime', 'idm', 'ibox'], plot_dir=plot_dir+str(gulp_i)+"_") # plot cluster results  