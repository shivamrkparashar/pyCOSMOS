#!/usr/bin/env python
import os
import operator
import sys
import string

path = r'Output/System_0'

out = open('half_out.dat','w')
out2 = open('full_out.dat','w')
#out.write("P [Pa]," + "P [atm]," + "A. Ads [cm^3/g]," + "Density [g/cm^3]," +"Total time [h],"+"Fug. Coeff.\n")
out.write("Pressure[atm]\t"  + "A.Ads[cm^3/g]\t"+"std\t"+ "Iso_Heat[kJ/mol]\t"+"std\n")
out2.write("Pressure[atm]\t" + "Pressure[kPa]\t" + "A.Ads[cm^3/g]\t"+"std\t"+ "A.Ads[molec/uc]\t"+"std\t"+
"Iso_Heat[kJ/mol]\t"+"std\t"+"Drift\t"+"TotalTime[h]\n")
results= []
results_full= []
#read the files in the folder "files"
for dir_entry in os.listdir(path):
    # Making all important variables = 1 to avoid errors when simulation is not finished.
    AAds_Molec_UC = str(1.0) ; AAds_Molec_UC_STD = str(1.0) 
    AAds_cm3_g    = str(1.0) ; AAds_cm3_g_STD= str(1.0)
    Drift         = str(1.0) ; TotalTimeSec  = str(1.0)
    IsoHeat       = str(1.0) ; IsoHeat_STD   = str(1.0)

    dir_entry_path = os.path.join(path, dir_entry)
    if os.path.isfile(dir_entry_path):
        with open(dir_entry_path, 'r') as my_file:
            for line in my_file:
                
                if "Partial pressure:" in line:
                    a,b,PartialPressure,c  = line.split()
                    #print PartialPressure[:-12]
                if "Fugacity coefficient" in line:
                    a,b,FugacityCoefficient,c  = line.split()
                    #print FucagityCoefficient
                #if "), density:" in line:
                    #a,b = line.split(", density:")
                    #print b
                    #c,d,AverageDensity,e = b.split()
                    #print AverageDensity[:-1]
                if "total time:" in line:
                    a,b,TotalTimeSec,c  = line.split()
                    #print TotalTimeSec
                if "Total energy-drift" in line:
                    a,b,Drift  = line.split()
                if "Enthalpy of adsorption:" in line:
                    for _ in range(9): my_file.next()
                    line = my_file.next()
#                    print line
#                    a,IsoHeat,b,IsoHeat_STD,c = line.split()
                    IsoHeat,b,IsoHeat_STD,c = line.split()
                if "Average loading absolute [molecules/unit cell]" in line:
                    a,b,c,d,e,AAds_Molec_UC,f,AAds_Molec_UC_STD,g = line.split()
                if "Average loading absolute [cm^3 (STP)/gr framework]" in line:
                    a,b,c,d,e,f,AAds_cm3_g,g,AAds_cm3_g_STD,h = line.split()

                    
#                    Heat of desorption:
#                    Total energy-drift:
        TotalTimeH = float(TotalTimeSec) / 3600
        #print TotalTimeH
        Pressure_atm = float(PartialPressure) / 101325
        Fugacity_bar = float(PartialPressure)*float(FugacityCoefficient)/100000
        #print Pressure_atm
        results.append((Pressure_atm,float(AAds_cm3_g),float(AAds_cm3_g_STD),-float(IsoHeat),float(IsoHeat_STD)))
        results_full.append((Pressure_atm,float(PartialPressure)/1000,float(AAds_cm3_g),float(AAds_cm3_g_STD),
                             float(AAds_Molec_UC),float(AAds_Molec_UC_STD),-float(IsoHeat),float(IsoHeat_STD),
                             float(Drift),TotalTimeH))

results.sort(key=lambda s: s[0], reverse=False)
for pair in results:
    out.write("{0[0]:<16.6e}{0[1]:> 12.6f}{0[2]:8.2f}{0[3]:> 16.5f}{0[4]:16.2f}\n".format(pair))
results_full.sort(key=lambda s: s[0], reverse=False)
for pair in results_full:
    out2.write("{0[0]:<16.4e}{0[1]:<16.4e}{0[2]:> 12.6f}{0[3]:8.2f}{0[4]:> 16.6f}{0[5]:8.2f}\
{0[6]:> 16.5f}{0[7]:16.2f}{0[8]:> 10.1e}{0[9]:10.2f}\n".format(pair))

out.close()
out2.close()
print 'Done'
#"C:/Program\ Files/Microsoft\ Office/Office15/EXCEL.exe" "out.csv"
#open out.csv -a C:/Program\ Files/Microsoft\ Office/Office15/EXCEL.exe
