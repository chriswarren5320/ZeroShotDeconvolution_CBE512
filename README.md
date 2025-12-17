# ZeroShotDeconvolution - characterizing model performance with various deconvolution regularization terms, benchmarked against gold standard deconvolution algorithms. 

This code is my implementation/modification of of a previously published "self-supervised deep-learning tool for instant denoising and super-resolution in optical fluorescence microscopy"

Credit is due to https://tristazeng.github.io/ZS-DeconvNet-page/, from which this code was modified.  

This package includes the Python implementation of training and inference, and MATLAB implementation of training data generation and simulation of PSF. 

This repository contains the code used to implement, modify, and test the ZS-Deconvnet denoising and deconvolution algorthim on simulated datasets and biological samples from our lab. 

 % ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



/data_augment_recorrupt_matlab</code> includes the MATLAB codes for generating training datasets and simulation of PSF 
/data_augment_recorrupt_matlab/GenData4ZS-DeconvNet` includes the MATLAB codes for generating training dataset for 2D and 3D ZS-DeconvNet, as well as one demo for generating simulated 3D PSF
/data_augment_recorrupt_matlab/XxUtils` includes common tool packages
/train_inference_python</code> includes the Python codes of training and inference, and the required dependencies
/train_inference_python/models</code> includes the optional models
/train_inference_python/utils</code> is the tool package

+ `./saved_models` includes pre-trained models for testing, and for each modality and structure `xx`:
  - `./saved_models/xx/saved_model` includes corresponding the pre-trained model and inference result, which you should be able to get by processing the raw test data with the pre-trained model
  - `./saved_models/xx/test_data` includes raw test data


My environment is:

- Windows 10
- Python 3.9.7
- Tensorflow-gpu 2.5.0
- NVIDIA GPU (GeForce RTX 3090) + CUDA (11.4)

To use our code, you should create a virtual environment and install the required packages first.

```
$ conda create -n zs-deconvnet python=3.9.7 
$ conda activate zs-deconvnet
$ pip install -r requirements.txt
```

After that, remember to install the right version of CUDA and cuDNN, if you want to use GPU. You can get the compatible version (e.g., cudatoolkit==11.3.1, cudnn==8.2.1) by searching

```
$ conda search cudatoolkit --info
$ conda search cudnn --info
```

then installing the corresponding version

```
$ conda install cudatoolkit==11.3.1
$ conda install cudnn==8.2.1
```


Dataset generation for ZS-DeconvNet:

+ Prepare a folder of raw data. Download [our open-source raw data with corresponding PSF](https://www.zenodo.org/record/7261163) of various modalities or use your own raw data. 

+ Open `./data_augment_recorrupt_matlab/GenData4ZS-DeconvNet/main_augm.m` and replace the parameter `data_folder` with your raw data directory. 

+ The default output path is `./your_augmented_datasets/`.



To train a new model, you need to:

+ Generated the training dataset following the instructions in the previous part and obtain the path to the corresponding PSF.

+ Choose a test image/volume.

+ Choose the demo file based on needs:
  
  + `./train_inference_python/train_demo_2D.sh`: train with 2D wide-field data and alike.
  


+ Set `otf_or_psf_path` (or `psf_path`), `data_dir`, `folder` and `test_images_path` in said demo file. Remember to examine other parameters like patch size and dx as well because their default value may not fit your training data.

+ Run it in your terminal.

+ The result wills be saved to <code>./your_saved_models/</code>.

+ Run <code>tensorboard --logdir [save_weights_dir]/[save_weights_name]/graph</code> to monitor the training process via tensorboard if needed.

+ Other **detailed description of each input argument of the python codes** can be found in the comments of the demo file. You can set them accordingly.

<hr>

<h2  id="Implementation of Python code1">5. Test a well-trained model</h2>

To test a well-trained ZS-DeconvNet model, you should:

+ Change the weight paths in <code>./train_inference_python/infer_demo_2D.sh</code> or <code>./train_inference_python/infer_demo_3D.sh</code> accordingly, or just use the default options given by us. The inference of SIM data is the same as other type of microscopy so no additional code is provided, though remember to set the model name correctly.
+ Run it in your terminal.
+ The output will be saved to the folder where you load weights, e.g., if you load weights from <code>./train_inference_python/saved_models/.../weights_40000.h5</code>, then the output will be saved to <code>./train_inference_python/saved_models/.../Inference/</code>.

<hr>
