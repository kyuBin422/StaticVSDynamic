clear;clc;close all
addpath(genpath('utils'));
format long
load matlab_x0=1point3.mat;

% calculate the expected omega and q
disp(mean(omegaStar(thetaqList,qList)))

disp(mean(qList))