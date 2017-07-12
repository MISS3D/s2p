import os
import numpy as np
import json
import pickle
import matplotlib.pyplot as plt
import sys
from PIL import Image
import subprocess

def read_paths_to_pointing_correction_txt_from_directory(main_directory):
	pointing_txt_path=[]
	for root,dirs,files in os.walk(main_directory):
		for file in files:
			if file.endswith("pointing.txt"):
				pointing_txt_path.append(os.path.join(root,file))
	return pointing_txt_path

def read_tile_pixelized_zone_of_imag_from_config_file_path(path_to_config_file):
	with open(path_to_config_file,'r') as f:
		cfg_local=json.load(f)
	min_row=(cfg_local['roi']['y'])/cfg_local['tile_size']
	max_row=(cfg_local['roi']['y']+cfg_local['roi']['h'])/cfg_local['tile_size']
	min_col=(cfg_local['roi']['x'])/cfg_local['tile_size']
	max_col=(cfg_local['roi']['x']+cfg_local['roi']['w'])/cfg_local['tile_size']
	return [min_row,max_row,min_col,max_col]

def get_translation_for_failling_tile(path_to_config_file,interpolated_translations):
	[min_row,max_row,min_col,max_col]=read_tile_pixelized_zone_of_imag_from_config_file_path(path_to_config_file)
	# [xtransltions,ytranslations]=get_interpolated_translationss(main)
	xtranslation=interpolated_translations[0][min_row][min_col]
	ytranslation=interpolated_translations[1][min_row][min_col]
	return np.array([[1, 0, xtranslation],[0, 1, ytranslation],[0, 0, 1]])

def read_tile_pixelized_zone_of_imag_from_config_file_of_main_directory(main_directory):
	return read_tile_pixelized_zone_of_imag_from_config_file_path(os.path.join(main_directory+'config.json'))
	# with open(os.path.join(main_directory+'config.json'),'r') as f:
	# 	cfg_local=json.load(f)
	# min_row=(cfg_local['roi']['y'])/cfg_local['tile_size']
	# max_row=(cfg_local['roi']['y']+cfg_local['roi']['h'])/cfg_local['tile_size']
	# min_col=(cfg_local['roi']['x'])/cfg_local['tile_size']
	# max_col=(cfg_local['roi']['x']+cfg_local['roi']['w'])/cfg_local['tile_size']
	# return [min_row,max_row,min_col,max_col]

#image_of_translations=np.zeros([125,45,2])
def produce_image_of_translations(main_directory,write_translations=False):
	pointing_txt_path=read_paths_to_pointing_correction_txt_from_directory(main_directory)
	image_of_translations=np.zeros([125,45,2])
	image_of_translations[:]=np.NAN
	for el in pointing_txt_path:
    		pathToJson=os.path.dirname(os.path.dirname(el))+"/config.json"
    		with open(pathToJson,'r') as f:
    		    a=json.load(f)
    		tile_size=a['tile_size']
    		col_position=int(a['roi']['x']/tile_size)
    		row_position=int(a['roi']['y']/tile_size)
    		#print(col_position);print(row_position)
    		matrix=np.loadtxt(el)
    		#print(matrix[0][2],matrix[1][2])
    		image_of_translations[row_position][col_position][0]=matrix[0][2]
    		image_of_translations[row_position][col_position][1]=matrix[1][2]
	if (write_translations==True):
		with open(os.path.join(main_directory,'translations.txt'),"wb") as f:
			pickle.dump(image_of_translations,f,protocol=2)
	return image_of_translations

def remove_NaN(main_directory,write_translations=False):
	"""remove NaN values by computing the approximate solutions of poisson or laplace equations for images"""
	with open(os.path.join(main_directory,'translations.txt'),"rb") as f:
			image_of_translations=pickle.load(f)
	xtranslation_path_tif=os.path.join(main_directory,"xtranslations.tif");xtranslation_no_nan_path_tif=os.path.join(main_directory,"xtranslations_no_nan.tif")
	ytranslation_path_tif=os.path.join(main_directory,"ytranslations.tif");ytranslation_no_nan_path_tif=os.path.join(main_directory,"ytranslations_no_nan.tif")
	img1=Image.fromarray(image_of_translations[:,:,0]);img1.save(xtranslation_path_tif)
	img2=Image.fromarray(image_of_translations[:,:,1]);img2.save(ytranslation_path_tif)
	subprocess.call(['../imscript/bin/simpois','-i',xtranslation_path_tif,'-o',xtranslation_no_nan_path_tif])
	subprocess.call(['../imscript/bin/simpois','-i',ytranslation_path_tif,'-o',ytranslation_no_nan_path_tif])

	images=[Image.open(xtranslation_no_nan_path_tif),Image.open(ytranslation_no_nan_path_tif)]
	image_of_translations_no_nan=[np.array(images[0],dtype=image_of_translations.dtype),np.array(images[1],dtype=image_of_translations.dtype)]
	if (write_translations==True):
		with open(os.path.join(main_directory,'translations_no_nan.txt'),"wb") as f:
			pickle.dump(image_of_translations_no_nan,f,protocol=2)
	return image_of_translations_no_nan

def compute_interpolated_translation(main_directory):
	produce_image_of_translations(main_directory,True)
	return remove_NaN(main_directory,True)

def get_interpolated_translations(main_directory):
	path_to_translations_no_nan=os.path.join(main_directory,'translations_no_nan.txt')
	with open(path_to_translations_no_nan,"rb") as f:
		image_of_translations_no_nan=pickle.load(f)	
	return image_of_translations_no_nan

def plot_translations_no_nan(main_directory):
	[min_row,max_row,min_col,max_col]=read_tile_pixelized_zone_of_imag_from_config_file_of_main_directory(main_directory)
	path_to_translations_no_nan=os.path.join(main_directory+'translations_no_nan.txt')
	with open(path_to_translations_no_nan,"rb") as f:
		image_of_translations_no_nan=pickle.load(f)
	plt.clf()
	plt.cla()
	plots=plt.subplots(1,2)
	plots[1][0].imshow(image_of_translations_no_nan[0][min_row:max_row,min_col:max_col]);plots[1][0].set_title("xtranslation with no NaN value")#+" for directory "+ main_directory)
	plots[1][1].imshow(image_of_translations_no_nan[1][min_row:max_row,min_col:max_col]);plots[1][1].set_title("ytranslation with no NaN value")#+" for directory "+ main_directory)
	#plots[1].set_title("x and y translation"+" for directory "+ main_directory)
	plt.show()

def plot_translations(main_directory):
	[min_row,max_row,min_col,max_col]=read_tile_pixelized_zone_of_imag_from_config_file_of_main_directory(main_directory)
	path_to_translations=os.path.join(main_directory+'translations.txt')
	with open(path_to_translations,"rb") as f:
		image_of_translations=pickle.load(f)
	plt.clf()
	plt.cla()
	plots=plt.subplots(1,2)
	plots[1][0].imshow(image_of_translations[min_row:max_row,min_col:max_col,0]);plots[1][0].set_title("xtranslation")#+" for directory "+ main_directory)
	plots[1][1].imshow(image_of_translations[min_row:max_row,min_col:max_col,1]);plots[1][1].set_title("ytranslation")#+" for directory "+ main_directory)
	#plots[1].set_title("x and y translation"+" for directory "+ main_directory)
	plt.show()

def plot_all_translations(main_directory):
	[min_row,max_row,min_col,max_col]=read_tile_pixelized_zone_of_imag_from_config_file_of_main_directory(main_directory)
	path_to_translations=os.path.join(main_directory+'translations.txt')
	path_to_translations_no_nan=os.path.join(main_directory+'translations_no_nan.txt')
	with open(path_to_translations,"rb") as f:
		image_of_translations=pickle.load(f)
	with open(path_to_translations_no_nan,"rb") as f:
		image_of_translations_no_nan=pickle.load(f)
	plt.clf()
	plt.cla()
	plots=plt.subplots(2,2)
	plots[1][0][0].imshow(image_of_translations[min_row:max_row,min_col:max_col,0]);plots[1][0][0].set_title("xtranslation")#+" for directory "+ main_directory)
	plots[1][0][1].imshow(image_of_translations[min_row:max_row,min_col:max_col,1]);plots[1][0][1].set_title("ytranslation")#+" for directory "+ main_directory)
	plots[1][1][0].imshow(image_of_translations_no_nan[0][min_row:max_row,min_col:max_col]);plots[1][1][0].set_title("xtranslation with no NaN value")#+" for directory "+ main_directory)
	plots[1][1][1].imshow(image_of_translations_no_nan[1][min_row:max_row,min_col:max_col]);plots[1][1][1].set_title("ytranslation with no NaN value")
	#plots[1].set_title("x and y translation"+" for directory "+ main_directory)
	plt.show()



if __name__=="__main__":
	main_directory=sys.argv[1]
	if (os.path.exists(os.path.join(main_directory,"translations.txt"))):
		print("the file "+os.path.join(main_directory,"translations.txt")+" already exists")
	else:
		print("computing the image of the pointing correction for each tile")	
		compute_interpolated_translation(main_directory)
	plot_all_translations(main_directory)
