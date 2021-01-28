%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Processamento de Sinal - BPM                          %
%                           13824  --  17618                              %
%                                                                         %
%   This script was made with the goal of get BPM's from 30 seconds of    %
%    some music applying fiters (if needed) to get the better result      %
%       possible and we script display the graphs that show the           %
%                    filter and the output signal                         %
%                                                                         %
%      You can use the script by calling the function get_BPMs that       %
% returns the output signal filtered, the sampling frequency and the bpm  %
%      The input arguments are the file (music) the interval and the      %   
%   the last give you the choice to show only the correlation plot or     %
%   multiple plots that gives you more information about the signal.      %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

% 1st music
fprintf('John Paesano On the Case: expected BPMs: 170\n');
[output_s1, Fs1, bpm1] = get_BPMs('john paesano_on the case.wav',0,0);

% 2nd music
fprintf('Nicky Jam, Will Smith, Era Istrefi_Live It Up : expected BPMs: 144\n');
[output_s2, Fs2, bpm2] = get_BPMs('nicky jam_live it up.wav',0,0);
 
% 3rd music
fprintf('Kodomo Concept 16 : expected BPMs: 94\n');
[output_s3, Fs3, bpm3] = get_BPMs('kodomo_concept 16.wav',0,0);
 
% 4th music
fprintf('Brenda Lee Rockin Around The Christmas Tree: expected BPMs: 67\n');
[output_s4, Fs4, bpm4] = get_BPMs('Brenda Lee_ rockin around.wav', [150 300],0);