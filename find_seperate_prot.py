from pypdb import *

def get_hit_params(hit):
    for data in hit:
        if "pdbid" in data: 
            pdb_name=data.split(":")[0]
            break
    for data in hit:
        if "Score" in data: break
    data=str(data)
    start=data.find("Score")
    end=data.find("Positives")
    eval= (data[start:end].split()[7][0:-1])
    identity=(data[start:end].split()[14][1:-3])
    return pdb_name,eval,identity


def findSeperate(pdb_id,chain_id):
    blast_results = get_blast2(pdb_id, chain_id=chain_id)
    i=0
    names=[]
    chain1=[]
    for i in range(len(blast_results[1])):
        if int(get_hit_params(blast_results[1][i])[2])>=95:
            names.append(get_hit_params(blast_results[1][i])[0])


    if len(names)==0:
        return "-1"
    for na in names:
        all_info = get_all_info(na)
        if len(all_info["polymer"])>7:
            chain1.append(na)
    if len(chain1)==0:
        return "-1"

    resolution = 100
    protein=chain1[0]
    if len(chain1)>1:
        for ch in chain1:
            result = describe_pdb(ch)
            if result['expMethod']=='X-RAY DIFFRACTION':
                if float(resolution)>float(result['resolution']):
                    resolution=result['resolution']
                    protein=ch
            else:
                return ch
    return protein


def seqSeperateProtein(pdb_id,chain_id,protein):
    blast_results = get_blast(pdb_id,chain_id)
    just_hits = blast_results['BlastOutput_iterations']['Iteration']['Iteration_hits']['Hit']
    i=0
    while just_hits[i]['Hit_def'][0:4]!=protein:
        i+=1
    return just_hits[i]['Hit_hsps']['Hsp']['Hsp_midline']

"""
protein1=findSeperate('1ACB','E')

for i in get_blast2('1A6W','L')[1]:
    print (get_hit_params(i))
#print (protein1)
"""
"""
all_info=get_all_info("12E8")
print(len(all_info["polymer"][0]['chain']))
list1=[ '1A37', '1A38', '1A3I', '1A3J', '1A3L', '1A3N', '1A3O', '1A4F', '1A4J', '1A4K', '1A4Y', '1A50', '1A5A', '1A5B', '1A5F', '1A5H', '1A5S', '1A6D', '1A6E', '1A6T', '1A6U', '1A6V', '1A6W', '1A6Z', '1A75', '1A7C', '1A7N', '1A7O', '1A7P', '1A7Q', '1A7R', '1A81', '1A9W', '1A9X', '1AA1', '1ABJ', '1ABO', '1ABR', '1ABW', '1ABY', '1ACB', '1ACM', '1AD0', '1AD9', '1AE6', '1AGR', '1AHG', '1AHJ', '1AI4', '1AI5', '1AI6', '1AI7', '1AIF', '1AIK', '1AIP', '1AJ7', '1AJ9', '1AJN', '1AJP', '1AJQ', '1AJS', '1AK4', '1AKS', '1ALL', '1ALX', '1ALZ', '1AM4', '1AN1', '1AOK', '1AON', '1AP2', '1APH', '1APM', '1APT', '1APU', '1APV', '1APW', '1APY', '1APZ', '1AQ7', '1AQC', '1AQK', '1ARO', '1AS4', '1AT1', '1ATN', '1ATP', '1AUI', '1AUS', '1AUT', '1AVA', '1AVF', '1AVO', '1AVP', '1AVW', '1AVX', '1AVZ', '1AW8', '1AWH', '1AWI', '1AWQ', '1AWR', '1AWS', '1AWT', '1AWU', '1AWV', '1AXC', '1AXD', '1AXI', '1AXS', '1AXT', '1AY1', '1AY7', '1AYA', '1AYB', '1AYC', '1AYY', '1AZZ', '1B05', '1B07', '1B0H', '1B0N', '1B17', '1B18', '1B19', '1B1H', '1B27', '1B2A', '1B2B', '1B2C', '1B2D', '1B2E', '1B2F', '1B2G', '1B2H', '1B2S', '1B2U', '1B2W', '1B32', '1B34', '1B3F', '1B3G', '1B3H', '1B3L', '1B3S', '1B40', '1B41', '1B46', '1B4H', '1B4J', '1B4U', '1B4Z', '1B51', '1B52', '1B58', '1B5F', '1B5H', '1B5I', '1B5J', '1B6C', '1B6H', '1B6J', '1B70', '1B7H', '1B7Y', '1B86', '1B8H', '1B8M', '1B9E', '1B9J', '1BAB', '1BAF', '1BBB', '1BBD', '1BBJ', '1BBZ', '1BC5', '1BCK', '1BDJ', '1BE9', '1BEN', '1BEU', '1BEY', '1BFO', '1BFV', '1BGS', '1BH8', '1BH9', '1BHF', '1BHH', '1BI4', '1BI7', '1BI8', '1BIJ', '1BIQ', '1BJ3', '1BJO', '1BJR', '1BKD', '1BKS', '1BLL', '1BLN', '1BLX', '1BM2', '1BM3', '1BMB', '1BML', '1BMQ', '1BND', '1BOU', '1BP3', '1BPH', '1BPL', '1BQM', '1BQN', '1BQP', '1BQQ', '1BR1', '1BR4', '1BR8', '1BRB', '1BRC', '1BRL', '1BRS', '1BS6', '1BS8', '1BSX', '1BT6', '1BUH', '1BUI', '1BUN', '1BUV', '1BUW', '1BVL', '1BVN', '1BVY', '1BW8', '1BW9', '1BWV', '1BX9', '1BXG', '1BXI', '1BXN', '1BXR', '1BXX', '1BY5', '1BZ0', '1BZ1', '1BZ7', '1BZH', '1BZQ', '1BZX', '1BZZ', '1C0F', '1C0G', '1C0T', '1C0U', '1C12', '1C16', '1C1B', '1C1C', '1C1D', '1C1E', '1C1X', '1C1Y', '1C29', '1C30', '1C3A', '1C3O', '1C40', '1C4Z', '1C5B', '1C5C', '1C5D', '1C5F', '1C5M', '1C5W', '1C5X', '1C5Y', '1C5Z', '1C6V', '1C7B', '1C7C', '1C7D', '1C8O', '1C8V', '1C9D', '1C9I', '1C9L', '1C9P', '1C9T', '1CA9', '1CAU', '1CAV', '1CAW', '1CAX', '1CB7', '1CC0', '1CC1', '1CCW', '1CD1', '1CD9', '1CDK', '1CDL', '1CDM', '1CE7', '1CE8', '1CF0', '1CF8', '1CFQ', '1CFV', '1CG5', '1CG8', '1CGI', '1CGJ', '1CGS', '1CI6', '1CIQ', '1CJF', '1CJQ', '1CJR', '1CK0', '1CKA', '1CKB', '1CKN', '1CLO', '1CLS', '1CLV', '1CLY', '1CLZ', '1CM1', '1CM4', '1CMI', '1CMK', '1CMX', '1CMY', '1CN3', '1CN4', '1CO7', '1COH', '1CP3', '1CP9', '1CPB', '1CPC', '1CPH', '1CPI', '1CQ4', '1CQI', '1CQJ', '1CR9', '1CS0', '1CSB', '1CSE', '1CSO', '1CT0', '1CT2', '1CT4', '1CT8', '1CTP', '1CVS', '1CVU', '1CVW', '1CW2', '1CWA', '1CWB', '1CWC', '1CWD', '1CWE', '1CWF', '1CWH', '1CWI', '1CWJ', '1CWK', '1CWL', '1CWM', '1CWO', '1CX9', '1CXP', '1CXZ', '1CYN', '1CZI', '1CZQ', '1CZY', '1CZZ', '1D00', '1D01', '1D09', '1D0A', '1D0D', '1D0G', '1D0J', '1D2V', '1D2Z', '1D3B', '1D4T', '1D4V', '1D4W', '1D4X', '1D5B', '1D5D', '1D5E', '1D5H', '1D5I', '1D5L', '1D5S', '1D6R', '1D6V', '1D6W', '1D7W', '1D8T', '1D9I', '1DBA', '1DBB', '1DBJ', '1DBK', '1DBM', '1DCE', '1DD3', '1DD4', '1DDV', '1DEI', '1DEJ', '1DEV', '1DF0', '1DFB', '1DFJ', '1DGH', '1DHK', '1DII', '1DIQ', '1DJ7', '1DJS', '1DKD', '1DKE', '1DKF', '1DKG', '1DKX', '1DKY', '1DKZ', '1DL7', '1DLF', '1DLO', '1DM0', '1DML', '1DN0', '1DN2', '1DNU', '1DNW', '1DOA', '1DOJ', '1DOW', '1DP5', '1DPH', '1DPJ', '1DQD', '1DQL', '1DQM', '1DQQ', '1DS2', '1DS5', '1DS6', '1DSF', '1DTD', '1DTQ', '1DTT', '1DTW', '1DU3', '1DXP', '1DXT', '1DXU', '1DXV', '1DY8', '1DY9', '1DZB', '1DZG', '1DZH', '1DZI', '1E0O', '1E1C', '1E1H', '1E3A', '1E3D', '1E3U', '1E3W', '1E44', '1E4E', '1E4K', '1E50', '1E54', '1E6E', '1E6I', '1E6O', '1E8N', '1E94', '1E96', '1E9H', '1E9O', '1E9P', '1E9Q', '1E9Y', '1E9Z', '1EAI', '1EAK', '1EAP', '1EAW', '1EAY', '1EBA', '1EBD', '1EBP', '1EE4', '1EE5', '1EEN', '1EEO', '1EER', '1EET', '1EF1', '1EFN', '1EFP', '1EFU', '1EFV', '1EG4', '1EG9', '1EGP', '1EJ4', '1EJ7', '1EJA', '1EJH', '1EJL', '1EJM', '1EJY', '1ELR', '1ELW', '1EM8', '1EMT', '1EMU', '1EMV', '1EO2', '1EO9', '1EOA', '1EOB', '1EOC', '1EOJ', '1EOL', '1EP1', '1EP2', '1EP3', '1EP4', '1EPL', '1EPM', '1EQY', '1ER8', '1ES0', '1ES7', '1ESV', '1ETH', '1ETR', '1ETS', '1ETT', '1ETZ', '1EUC', '1EUD', '1EUI', '1EUV', '1EV2', '1EV3', '1EV6', '1EVH', '1EVR', '1EVT', '1EWY', '1EXB', '1EXU', '1EZQ', '1EZS', '1EZU', '1EZZ', '1F02', '1F0C', '1F0R', '1F0S', '1F11', '1F1B', '1F1J', '1F1W', '1F2S', '1F2T', '1F2U', '1F34', '1F3D', '1F3M', '1F3U', '1F3V', '1F45', '1F47', '1F4V', '1F4W', '1F4X', '1F4Y', '1F51', '1F59', '1F5Q', '1F5R', '1F60', '1F6A', '1F6F', '1F6M', '1F7A', '1F7Z', '1F80', '1F8A', '1F8T', '1F8U', '1F93', '1F99', '1FAI', '1FAP', '1FAV', '1FAW', '1FAX', '1FC2', '1FCC', '1FCD', '1FCH', '1FDH', '1FEV', '1FFG', '1FFS', '1FFW', '1FG9', '1FGL', '1FGN', '1FGV', '1FH5', '1FHJ', '1FIG', '1FIN', '1FIP', '1FIV', '1FIW', '1FIZ', '1FJM', '1FJS', '1FK9', '1FKN', '1FKO', '1FKP', '1FL3', '1FL5', '1FL6', '1FL7', '1FLC', '1FLE', '1FLL', '1FLR', '1FLT', '1FM0', '1FM2', '1FMA', '1FMO', '1FN3', '1FN4', '1FN8', '1FNE', '1FNG', '1FOE', '1FOR', '1FP4', '1FPP', '1FPR', '1FQ1', '1FQ9', '1FQK', '1FQV', '1FR2', '1FRF', '1FRV', '1FS0', '1FS1', '1FS2', '1FSS', '1FSX', '1FT1', '1FT2', '1FUY', '1FVC', '1FVD', '1FVE', '1FVM', '1FVU', '1FVV', '1FX0', '1FXH', '1FXV', '1FXW', '1FY4', '1FY5', '1FY8', '1FYH', '1FYN', '1FYR', '1G08', '1G09', '1G0A', '1G0B', '1G0V', '1G0Y', '1G1F', '1G1G', '1G1H', '1G1S', '1G2C', '1G2L', '1G2M', '1G37', '1G3I', '1G3J', '1G4A', '1G4B', '1G4U', '1G4Y', '1G5Q', '1G6G', '1G6V', '1G72', '1G73', '1G7A', '1G7B', '1G7C', '1G8J', '1G8K', '1G9I', '1G9V', '1G9W', '1GA1', '1GA4', '1GA6', '1GAF', '1GAG', '1GAQ', '1GBB', '1GBC', '1GBD', '1GBF', '1GBH', '1GBI', '1GBK', '1GBL', '1GBM', '1GBU', '1GBV', '1GCQ', '1GCV', '1GCW', '1GDN', '1GDQ', '1GDU', '1GEC', '1GFW', '1GGB', '1GGC', '1GGP', '1GH0', '1GH6', '1GHD', '1GHF', '1GHQ', '1GI7', '1GI8', '1GI9', '1GIG', '1GJ7', '1GJ8', '1GJ9', '1GJA', '1GJB', '1GJC', '1GJD', '1GK0', '1GK1', '1GK8', '1GK9', '1GKA', '1GKF', '1GL0', '1GL1', '1GL4', '1GLA', '1GLB', '1GLC', '1GLD', '1GLE', '1GLI', '1GM7', '1GM8', '1GM9', '1GMW', '1GNG', '1GO3', '1GO4', '1GO6', '1GPO', '1GPQ', '1GPW', '1GRN', '1GTJ', '1GTL', '1GUA', '1GUJ', '1GUL', '1GVK', '1GVN', '1GVU', '1GVX', '1GWQ', '1GWR', '1GXC', '1GXD', '1GXS', '1GYB', '1GZL', '1GZP', '1GZQ', '1GZS', '1GZW', '1GZX', '1H0G', '1H0H', '1H0I', '1H1L', '1H1P', '1H1Q', '1H1R', '1H1S', '1H1V', '1H2A', '1H2G', '1H2K', '1H2L', '1H2M', '1H2R', '1H2S', '1H2T', '1H2U', '1H2V', '1H31', '1H32', '1H33', '1H3O', '1H3P', '1H4I', '1H4J', '1H4L', '1H59', '1H5R', '1H5T', '1H6E', '1H6K', '1H6W', '1H9H', '1H9I', '1H9L', '1H9O', '1H9S', '1HA7', '1HAB', '1HAC', '1HAG', '1HAI', '1HAV', '1HAX', '1HAZ', '1HBA', '1HBB', '1HBH', '1HBR', '1HBS', '1HCF', '1HCG', '1HCN', '1HCO', '1HDA', '1HDB', '1HDM', '1HDS', '1HE1', '1HE8', '1HEF', '1HES', '1HFE', '1HGA', '1HGB', '1HGC', '1HGD', '1HGE', '1HGF', '1HGG', '1HGH', '1HGI', '1HGJ', '1HH4', '1HHO', '1HHZ', '1HIL', '1HKD', '1HKL', '1HL3', '1HL6', '1HLA', '1HLE', '1HLU', '1HMV', '1HNE', '1HNI', '1HNV', '1HPG', '1HPZ', '1HQ4', '1HQ6', '1HQE', '1HQQ', '1HQU', '1HQW', '1HQY', '1HR6', '1HR7', '1HRP', '1HT1', '1HT2', '1HTE', '1HTM', '1HTR', '1HTV', '1HUC', '1HV4', '1HWG', '1HWH', '1HWM', '1HWN', '1HWO', '1HWP', '1HX1', '1HXL', '1HXM', '1HXZ', '1HY2', '1HYR', '1HZH', '1HZZ', '1I1Q', '1I1R', '1I2M', '1I31', '1I3G', '1I3R', '1I3Z', '1I4D', '1I4E', '1I4L', '1I4O', '1I4T', '1I5K', '1I5O', '1I72', '1I73', '1I79', '1I7A', '1I7B', '1I7C', '1I7M', '1I7Q', '1I7S', '1I7W', '1I7X', '1I7Y', '1I7Z', '1I85', '1I8L', '1I9C', '1I9I', '1I9J', '1IAO', '1IAR', '1IAU', '1IB1', '1IBE', '1IBG', '1IBR', '1IBT', '1IBU', '1IBV', '1IBW', '1IEA', '1IEB', '1IGF', '1IGI', '1IGJ', '1IGM', '1IGT', '1IGY', '1IHJ', '1II4', '1II8', '1IID', '1IIL', '1IJD', '1IJE', '1IJF', '1IK9', '1IKV', '1IKW', '1IKX', '1IKY', '1IL1', '1IND', '1INE', '1IOE', '1IQ1', '1IQ5', '1IQE', '1IQF', '1IQG', '1IQH', '1IQI', '1IQJ', '1IQK', '1IQL', '1IQM', '1IQN', '1IQW', '1IR1', '1IR2', '1IR3', '1IRA', '1IRD', '1IRE', '1IS0', '1IS7', '1IS8', '1ISQ', '1IT9', '1ITB', '1IVO', '1IWA', '1IWH', '1IWQ', '1IXR', '1IXS', '1IXX', '1IYJ', '1IZA', '1IZB', '1IZE', '1IZN', '1J05', '1J19', '1J2J', '1J2Q', '1J2X', '1J3I', '1J3J', '1J3K', '1J3Y', '1J3Z', '1J40', '1J41', '1J4X', '1J71', '1J73', '1J7D', '1J7S', '1J7V', '1J7W', '1J7Y', '1J7Z', '1J80', '1J81', '1J82', '1JAC', '1JAN', '1JAP', '1JAT', '1JBO', '1JBP', '1JCA', '1JCH', '1JCK', '1JCQ', '1JD5', '1JD6', '1JDB', '1JDH', '1JDP', '1JEB', '1JEK', '1JEN', '1JEQ', '1JET', '1JEU', '1JEV', '1JFQ', '1JG3', '1JGL', '1JGU', '1JGV', '1JHK', '1JIW', '1JJC', '1JK0', '1JK4', '1JK9', '1JKG', '1JKH', '1JKJ', '1JKY', '1JLA', '1JLB', '1JLC', '1JLE', '1JLF', '1JLG', '1JLL', '1JLQ', '1JLT', '1JLU', '1JMA', '1JMT', '1JN6', '1JN9', '1JNH', '1JNL', '1JNN', '1JNR', '1JNZ', '1JOJ', '1JOT', '1JOU', '1JOW', '1JP5', '1JPL', '1JPP', '1JPT', '1JPW', '1JQ8', '1JQ9', '1JQJ', '1JQL', '1JRO', '1JRP', '1JRR', '1JRS', '1JRT', '1JSD', '1JSH', '1JSI', '1JSM', '1JSN', '1JSO', '1JST', '1JTD', '1JTG', '1JTH', '1JTO', '1JTP', '1JTT', '1JUI', '1JUQ', '1JV2', '1JV5', '1JVZ', '1JW0', '1JW6', '1JW9', '1JWA', '1JWB', '1JWG', '1JWH', '1JWI', '1JX9', '1JXP', '1JXQ', '1JY7', '1JYC', '1JYI', '1JYO', '1JYQ', '1JYR', '1JZD', '1K0Y', '1K1K', '1K28', '1K2X', '1K3A', '1K3U', '1K3Z', '1K4W', '1K54', '1K55', '1K56', '1K57', '1K5Q', '1K5S', '1K6Q', '1K7D', '1K7E', '1K7F', '1K7L', '1K7X', '1K8I', '1K8R', '1K8X', '1K8Y', '1K8Z', '1K90', '1K93', '1K9O', '1KA9', '1KAC', '1KAP', '1KC2', '1KCG', '1KCU', '1KCV', '1KD2', '1KD8', '1KD9', '1KDD', '1KDQ', '1KDV', '1KDY', '1KDZ', '1KE1', '1KE2', '1KEC', '1KEE', '1KEL', '1KEM', '1KF9', '1KFA', '1KFB', '1KFC', '1KFE', '1KFJ', '1KFK', '1KFU', '1KFX', '1KGC', '1KGY', '1KHP', '1KHQ', '1KI1', '1KI6', '1KIU', '1KJ1', '1KJ4', '1KJ7', '1KJF', '1KJG', '1KJH', '1KJY', '1KKL', '1KKM', '1KKQ', '1KL3', '1KL5', '1KLF', '1KLI', '1KLJ', '1KLM', '1KMC', '1KMH', '1KMI', '1KN1', '1KN2', '1KN4', '1KNA', '1KNE', '1KNO', '1KO6', '1KPS', '1KRL', '1KSG', '1KSH', '1KSJ', '1KSN', '1KT2', '1KTD', '1KTK', '1KTP', '1KTZ', '1KU6', '1KU8', '1KUG', '1KUI', '1KUJ', '1KUK', '1KV6', '1KVD', '1KVE', '1KXP', '1KXQ', '1KXT', '1KXV', '1KY6', '1KY7', '1KYD', '1KYF', '1KYI', '1KYU', '1KZ7', '1KZG', '1KZU', '1KZY', '1KZZ', '1L0A', '1L0O', '1L0X', '1L0Y', '1L2I', '1L2W', '1L3R', '1L4D', '1L4Z', '1L5H', '1L6O', '1L6X', '1L7I', '1L7T', '1L7V', '1L7Z', '1LA6', '1LB1', '1LB5', '1LB6', '1LCJ', '1LCK', '1LD7', '1LD8', '1LDJ', '1LDT', '1LEM', '1LEN', '1LES', '1LEW', '1LEZ', '1LF8', '1LFD', '1LFL', '1LFQ', '1LFT', '1LFV', '1LFY', '1LFZ', '1LGH', '1LI1', '1LIA', '1LJ2', '1LJW', '1LKK', '1LKL', '1LKY', '1LMW', '1LNU', '1LO0', '1LO2', '1LO3', '1LO4', '1LOA', '1LOB', '1LOD', '1LOE', '1LOG', '1LOP', '1LOT', '1LP1', '1LPA', '1LPB', '1LPG', '1LPH', '1LPK', '1LPZ', '1LQ8', '1LQD', '1LQG', '1LQM', '1LQS', '1LQV', '1LRW', '1LS5', '1LSH', '1LT3', '1LT4', '1LUC', '1LUJ', '1LVB', '1LVC', '1LW0', '1LW2', '1LW6', '1LWC', '1LWE', '1LWF', '1LX5', '1LYA', '1LYW', '1LZW', '1M10', '1M1D', '1M1E', '1M1N', '1M1X', '1M21', '1M26', '1M2O', '1M2T', '1M2V', '1M2Z', '1M3D', '1M43', '1M45', '1M46', '1M4H', '1M4U', '1M5A', '1M5N', '1M6V', '1M71', '1M72', '1M7D', '1M7E', '1M7I', '1M9C', '1M9D', '1M9E', '1M9F', '1M9P', '1M9X', '1M9Y', '1MA3', '1MA9', '1MAE', '1MAF', '1MAH', '1MAM', '1MBU', '1MBV', '1MBX', '1MCB', '1MCC', '1MCD', '1MCE', '1MCF', '1MCH', '1MCI', '1MCJ', '1MCK', '1MCL', '1MCN', '1MCO', '1MCP', '1MCQ', '1MCR', '1MCS', '1MCT', '1MCV', '1MCW', '1MDU', '1MEE', '1MEL', '1MEX', '1MF2', '1MF4', '1MFA', '1MFB', '1MFC', '1MFD', '1MFE', '1MFG', '1MFL', '1MG9', '1MH2', '1MH5', '1MHL', '1MHM', '1MIE', '1MIK', '1MIM', '1MIO', '1MIU', '1MIZ', '1MJ7', '1MJ8', '1MJG', '1MJJ', '1MJU', '1MK2', '1MK7', '1MK9', '1MKO', '1ML0', '1MLB', '1MNF', '1MNU', '1MO1', '1MOX', '1MPJ', '1MQ5', '1MQ6', '1MQ8', '1MQK', '1MQL', '1MQM', '1MQN', '1MQS', '1MR1', '1MRC', '1MRD', '1MRE', '1MRF', '1MSO', '1MT1', '1MT7', '1MT8', '1MT9', '1MTP', '1MU2', '1MV9', '1MVC', '1MVF', '1MXE', '1MYP', '1MZ8', '1MZC', '1MZN', '1MZW']
print (len(list1))
complex2file=open("complex2file.txt","w")

for l in list1:
    complex2file.write(l)
    complex2file.write("\n")
"""
