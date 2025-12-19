
# ------------------------------- inference arguments -------------------------------

# --bs: batch size, a tuple of length n for n input images. e.g. 3 7 for 2 input images
# --num_seg_window_x: number of segmented patch along x axis for each input image, a tuple of length n for n input images. e.g. 3 7 for 2 input images
# --num_seg_window_y: number of segmented patch along y axis for each input image, a tuple of length n for n input images
# --overlap_x: overlapping in x direction for each input image, a tuple of length n for n input images. NOTICE: if num_seg_window>1, overlap should be big enough to avoid segmentation gap
# --overlap_y: overlapping in y direction for each input image, a tuple of length n for n input images.
# --predict_iter: loading weights iterations
# --input_dir: the path of test data
# --load_weights_path: root directory where models weights are loaded
# --insert_xy: padded blank margin in pixels. The applied blank margin may be larger than your given value due to size requirement of network. 
# --upsample_flag: 0 or 1, whether the network upsamples the image

# ------------------------------- examples -------------------------------

# WF lysosome
python Infer_2D.py --input_dir 'E:\Christopher\ZSD\train_inference_python\AAA_Cells_Selected\*.tif' \
                   --load_weights_path "E:\Christopher\ZSD\train_inference_python\20251119\weights\tubules_62_nosample_bg_700_std_100_SegNum2000_twostage_Unet_dxypsf_62_upsample_0_Hess_0.02\weights_10000.h5" \
                   --upsample_flag 0

