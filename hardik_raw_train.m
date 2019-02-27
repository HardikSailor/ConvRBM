clear all; clc; close all;
%%% Auditory Filterbank Learning Using ConvRBM.
%%% The code is developed by Hardik B. Sailor from ConvRBM code developed
%%% by H. Lee for spectrograms. We have developed it to learn auditory
%%% filterbank from the variable length raw speech and audio signals. We
%%% also used Noisy ReLU for inference, annealing dropout and Adam optimization
%%% First parameter is ws = 176 samples (8 ms) for sampling frequency =
%%% 22kHz (8 ms X 22 kHz)
%%% Second parameter is num_bases = 60. Number of filters
%%% epsilon = 0.001 learning rate paramneter, when you observe the weights
%%% diverge then set it to lower value such as 0.0005.
%%% l2reg = 0.001 weight decay parameter
%%% epsdecay = 0.0001 only used with SGD training. Currently using Adam
test2_ReLU_RAW_CRBM_hardik_main(176, 60, 0.001, 0.001, 0.0001);

%%% After ConvRBM training learned weights can be visualized as follows
%%% multiple_subplots(10,6,W,60) where W is ConvRBM weights
%%% You can sort weights as given in script filter_chars.m and see using
%%% above script
%%% After training ConvRBM, learned weights (filterbank) can be used to
%%% extract short-time features as an example given in TIMIT_raw_speech.m

