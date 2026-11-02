% The script run all validation tests implemented in the class Validate

clear;
clc;
close all;

% Add src to path
path(path, ['..', filesep, 'src']);

% Run all tests
[test] = Validate.allTests();