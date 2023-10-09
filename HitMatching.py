import numpy as np
import matplotlib as plt 
import sys, getopt


#compare absRS hits to standard hits to get some form of hit properties mapping between the two. 

def main(argv):
    inabs = ''
    inraw_p = ''
    inraw_n = ''  
    outfile = ''
    try:
        opts, args = getopt.getopt(argv, "ha:p:n:o:", ["help=","abs=","positive=","negative=","output="])
    except getopt.GetoptError:
        print("HitMatching.py -a <abs data in> -p <positive raw data in>  -n <negative raw data in> -o <outputfile>")
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print("HitMatching.py -a <abs data in> -p <positive raw data in>  -n <negative raw data in> -o <outputfile>")
            sys.exit()
        elif opt in ("-a", "--abs"):
            inabs = arg
        elif opt in ("-p", "--positive"):
            inraw_p = arg
        elif opt in ("-p", "--negative"):
            inraw_n = arg
        elif opt in ("-o", "--output"):
            outfile = arg

    #Hit data
    raw_data_p = np.loadtxt(inraw_p)
    raw_data_n = np.loadtxt(inraw_n)
    abs_data = np.loadtxt(inabs)

    filename = outfile


    class Hit:
        HitStart = -999
        HitEnd = -999
        TOT = -999
        PeakT = -999
        PeakADC = -999
        SADC = -999

        n_events = 501

    with open(filename, 'w') as f_out:
        
        #Run over events 
        for e in range(1,n_events):
            hit_window_sig = 5

            #Run over abs hits 
            for i in range(0,len(abs_data)):

                #abs hits for current event
                if (abs_data[i][0] == e):

                    #current abs hit
                    e, ch, startT, endT, peakT, TOT, PeakADC, SADC = abs_data[i]
                    n_abshits +=1

                    #'standard' hits for the +ve and -ve parts of the waveform
                    hit_neg = Hit()
                    hit_pos = Hit() 

                    #run over the original +ve hit data
                    for j in range(0,len(raw_data)):

                        #+ve raw hits for current event
                        if (raw_data[j][0] == e):

                            #current hit
                            m_e, m_ch, m_startT, m_endT, m_peakT, m_TOT, m_PeakADC, m_SADC = raw_data[j]

                            #check if channels match
                            if (m_ch == ch):

                                #see if +ve peak nested within abs hit:
                                if (m_startT > startT-hit_window_sig and m_endT < endT + hit_window_sig):

                                    hit_pos.HitStart = m_startT
                                    hit_pos.HitEnd = m_endT
                                    hit_pos.TOT = m_TOT
                                    hit_pos.PeakT = m_peakT
                                    hit_pos.PeakADC = m_PeakADC
                                    hit_pos.SADC = m_SADC

                    #run over the original -ve hit data
                    for k in range(0,len(raw_data_neg)):
                                        
                        #-ve raw hits for current event 
                        if (raw_data_neg[k][0] == e):

                            #current hit
                            m_e, m_ch, m_startT, m_endT, m_peakT, m_TOT, m_PeakADC, m_SADC = raw_data_neg[k]

                            #check if channels match
                            if (m_ch == ch):

                    
                                #see if -ve peak nested within abs hit:
                                if (m_startT > startT-hit_window_sig and m_endT < endT + hit_window_sig):

                                    hit_neg.HitStart = m_startT
                                    hit_neg.HitEnd = m_endT
                                    hit_neg.TOT = m_TOT
                                    hit_neg.PeakT = m_peakT
                                    hit_neg.PeakADC = m_PeakADC
                                    hit_neg.SADC = m_SADC
                                    
                    # if both, +ve and -ve peaks matched to abs hit: save the information for comparison
                    if (hit_neg.PeakT != -999 and hit_pos.PeakT != -999):
                        n_matched_abshits += 1
                        print_info = False 
                        if print_info:
                            print("---------------------------------------------------------------------------------")
                            print( "abs hit: \n", '{} {} {} {} {} {}'.format(ch, startT, endT, TOT, PeakADC, SADC))
                            print( "+ve hit: \n", '{} {} {} {} {} {}'.format(ch, hit_pos.HitStart, hit_pos.HitEnd,
                                                                             hit_pos.TOT, hit_pos.PeakT, 
                                                                             hit_pos.PeakADC, hit_pos.SADC))
                            print( "-ve hit: \n", '{} {} {} {} {} {}'.format(ch, hit_neg.HitStart, hit_neg.HitEnd,
                                                                             hit_negeg.TOT, hit_neg.PeakT, 
                                                                             hit_neg.PeakADC, hit_neg.SADC))
                        f_out.write('{} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {} {}\n'.format(e,ch, startT, endT, TOT, peakT,
 PeakADC, SADC,
                                                                    hit_pos.HitStart, hit_pos.HitEnd, hit_pos.TOT, hit_pos.PeakT, hit_p
os.PeakADC, hit_pos.SADC,
                                                                    hit_neg.HitStart, hit_neg.HitEnd, hit_neg.TOT, hit_neg.PeakT, hit_n
eg.PeakADC, hit_neg.SADC))
                        
    print(n_abshits, n_matched_abshits, n_matched_abshits/n_abshits )



if __name__ == "__main__":
    main(sys.argv[1:])

