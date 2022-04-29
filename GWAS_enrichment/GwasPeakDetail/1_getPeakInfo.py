#!/home/ajl1213/anaconda2/bin/python


import os
import operator



DIR=os.getcwd()
inFile='/home/ajl1213/Projects/Endothelial/Analysis2/GwasEnrich/1_getGwasCatalog/RESULT/MergedEC.SigGwasList.txt'
peakIntersectDIR='/home/ajl1213/Projects/Endothelial/Analysis2/GwasEnrich/1_getGwasCatalog/ObservedData'
diffPeakDIR='/home/ajl1213/Projects/Endothelial/data/ChIP/Processing/PeakBED'
cellTypeList=['EarlyEC','MidEC','LateEC','FullEC']

curateDict={
# Disease
'Coronary_artery_disease':['Disease','Yes'],
'Multiple_sclerosis':['Disease','Yes'],
'Asthma':['Disease','Yes'],
'Inflammatory_bowel_disease':['Disease','Yes'],
'Colorectal_cancer':['Disease','Yes'],
'Atrial_fibrillation':['Disease','Yes'],
'Crohns_disease':['Disease','Yes'],
'Eczema':['Disease','Yes'],
'Rheumatoid_arthritis':['Disease','Yes'],
'Ulcerative_colitis':['Disease','Yes'],
'Chronic_inflammatory_diseases_ankylosing_spondylitis,_Crohns_disease,_psoriasis,_primary_sclerosing_cholangitis,_ulcerative_colitis_pleiotropy':['Disease','No'],
'Autism_spectrum_disorder_or_schizophrenia':['Disease','No'],
'Atopic_asthma':['Disease','Yes'],
'Allergic_disease_asthma,_hay_fever_or_eczema':['Disease','No'],
'Alzheimers_disease_late_onset':['Disease','Yes'],
'Varicose_veins':['Disease','Yes'],

# Hematological trait
'Platelet_count':['HematologicalTrait','Yes'],
'Red_blood_cell_count':['HematologicalTrait','Yes'],
'White_blood_cell_count':['HematologicalTrait','Yes'],
'Lymphocyte_counts':['HematologicalTrait','Yes'],
'Red_cell_distribution_width':['HematologicalTrait','No'],
'Mean_platelet_volume':['HematologicalTrait','No'],
'Monocyte_count':['HematologicalTrait','Yes'],
'Neutrophil_count':['HematologicalTrait','Yes'],
'Hematocrit':['HematologicalTrait','Yes'],
'Platelet_distribution_width':['HematologicalTrait','No'],
'Mean_spheric_corpuscular_volume':['HematologicalTrait','No'],
'Plateletcrit':['HematologicalTrait','Yes'],
'Monocyte_percentage_of_white_cells':['HematologicalTrait','No'],
'Reticulocyte_fraction_of_red_cells':['HematologicalTrait','No'],
'Neutrophil_percentage_of_white_cells':['HematologicalTrait','No'],
'Eosinophil_percentage_of_white_cells':['HematologicalTrait','No'],
'Lymphocyte_percentage_of_white_cells':['HematologicalTrait','No'],
'Mean_corpuscular_hemoglobin_concentration':['HematologicalTrait','No'],
'Basophil_count':['HematologicalTrait','No'],
'Immature_fraction_of_reticulocytes':['HematologicalTrait','No'],
'Platelet-to-lymphocyte_ratio':['HematologicalTrait','No'],
'Basophil_percentage_of_white_cells':['HematologicalTrait','No'],

# Other trait
'Heel_bone_mineral_density':['OtherTrait','No'],
'Systolic_blood_pressure':['OtherTrait','Yes'],
'Waist-to-hip_ratio_adjusted_for_BMI':['OtherTrait','No'],
'Pulse_pressure':['OtherTrait','Yes'],
'Lung_function_FEV1_FVC':['OtherTrait','No'],
'Waist-hip_ratio':['OtherTrait','No'],
'Serum_total_protein_level':['OtherTrait','No'],
'Hip_index':['OtherTrait','No'],
'Total_cholesterol_levels':['OtherTrait','Yes'],
'Aspartate_aminotransferase_levels':['OtherTrait','No'],
'Waist-hip_index':['OtherTrait','No'],
'Apolipoprotein_A1_levels':['OtherTrait','Yes'],
'Electrocardiogram_morphology_amplitude_at_temporal_datapoints':['OtherTrait','No'],
'Triglycerides':['OtherTrait','Yes'],
'PR_interval':['OtherTrait','Yes'],
'Alanine_aminotransferase_levels':['OtherTrait','No'],
'Calcium_levels':['OtherTrait','No'],
'A_body_shape_index':['OtherTrait','No'],
'Gamma_glutamyl_transpeptidase':['OtherTrait','No'],
'Birth_weight':['OtherTrait','No'],
'Mean_arterial_pressure':['OtherTrait','Yes'],
'LDL_cholesterol':['OtherTrait','Yes'],
'Blood_urea_nitrogen_levels':['OtherTrait','No'],
'Bioavailable_testosterone_levels':['OtherTrait','No'],
'Glaucoma_primary_open-angle':['OtherTrait','No'],
'Serum_albumin_levels':['OtherTrait','No'],
'Gamma_glutamyl_transferase_levels':['OtherTrait','No'],
'Non-albumin_protein_levels':['OtherTrait','No'],
'Blond_vs._brown_black_hair_color':['OtherTrait','No'],
'LDL_cholesterol_levels':['OtherTrait','No'],
'Liver_enzyme_levels_gamma-glutamyl_transferase':['OtherTrait','No'],
'Fasting_glucose':['OtherTrait','No'],
'Apolipoprotein_B_levels':['OtherTrait','Yes'],
'Testosterone_levels':['OtherTrait','No'],
'Snoring':['OtherTrait','No'],
'Medication_use_vasodilators_used_in_cardiac_diseases':['OtherTrait','No'],
'Blood_glucose_levels':['OtherTrait','No'],
'Total_triglycerides_levels':['OtherTrait','Yes']
}

os.system('mkdir '+DIR+'/DataLabel')
os.system('mkdir '+DIR+'/GwasPeakData')



def attachCategory():
    input1=open(inFile,'r')
    all_input1=input1.readlines()
    output1=open(DIR+'/DataLabel/MergedEC.SigGwasList.Label.txt','w')
    output1.write('gwasID\tcategoryID\tfigLabel\tnPeak\tfcVal\tempiricalPval\tadjPval\n')
    for line in all_input1[1:]:
	[gwasID, nPeak, fcVal, empiricalPval, adjPval]=line.strip().split('\t')

	categoryID=curateDict[gwasID][0]
	figLabel=curateDict[gwasID][1]

	newLine=[gwasID, categoryID, figLabel, nPeak, fcVal, empiricalPval, adjPval]
	output1.write('\t'.join(newLine)+'\n')

    input1.close()


def loadGwasPeak():
    gwasPeakDict={}
    nPeakDict={}
    for gwasID in curateDict:
	figLabel=curateDict[gwasID][1]
	if figLabel == 'Yes':
	    input1=open(peakIntersectDIR+'/'+gwasID+'.MergedEC.Peak.intersect.bed', 'r')
	    all_input1=input1.readlines()
	    nPeak=0
	    for line in all_input1:
		each=line.strip().split('\t')
		chrID=each[0]
		pt1=each[1]
		pt2=each[2]
		peakID=chrID+':'+pt1+'-'+pt2
		nPeak+=1
		if gwasPeakDict.has_key(gwasID):
		    tmp=gwasPeakDict[gwasID]
		    tmp.append(peakID)
		    gwasPeakDict[gwasID]=tmp
		else:
		    gwasPeakDict[gwasID]=[peakID]
	    nPeakDict[gwasID]=nPeak
	    input1.close()

    return nPeakDict, gwasPeakDict


def loadDiffPeak():
    peakDict={}
    for cellType in cellTypeList:
	input1=open(diffPeakDIR+'/'+cellType+'.bed','r')
	all_input1=input1.readlines()
	for line in all_input1:
	    each=line.strip().split('\t')
	    peakID=each[3]
	    peakDict[peakID]=cellType
	input1.close()
	
    return peakDict


def makeData():
    nPeakDict, gwasPeakDict = loadGwasPeak()
    peakDict = loadDiffPeak()

    output1=open(DIR+'/GwasPeakData/GwasPeakData.txt','w')
    output1.write('traitID\tgwasID\tcellType\tnPeakTotal\tnPeak\n')
    for i in ['Disease','HematologicalTrait','OtherTrait']:
	for key in sorted(nPeakDict.items(), key=operator.itemgetter(1), reverse=True):
	    gwasID=key[0]
	    nPeakTotal=key[1]
	    traitID=curateDict[gwasID][0]
	    if traitID == i:

		tmpDict={}
		for peakID in gwasPeakDict[gwasID]:
		    cellType=peakDict[peakID]
		    if tmpDict.has_key(cellType):
			tmpDict[cellType]+=1
		    else:
			tmpDict[cellType]=1

		for cellType in cellTypeList:
		    if tmpDict.has_key(cellType):
			nPeak=tmpDict[cellType]
		    else:
			nPeak=0

		    newLine=[i, gwasID, cellType, str(nPeakTotal), str(nPeak)]
		    output1.write('\t'.join(newLine)+'\n')
    output1.close()




attachCategory()
makeData()



