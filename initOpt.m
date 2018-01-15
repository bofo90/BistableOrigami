function opt=initOpt(varargin)
%initOpt() initialize the simulations with default options. To change 
%specific options use initOpt('option1',value1,'option2',value2,...).
%
%The most important options are:
%1) 'inputType': 'individual'
%2) 'template': any string indicating the specific template used to generate 
%   the architected material. 
%   -> When 'individual' is used, possible values are the names of the
%   platonic solids, archimedean solids and some prisms, including: 
%   'tetrahedron', 'cube', 'octahedron', truncated tetrahedron',
%   'cuboctahedron', 'truncated cube', 'truncated octahedron',
%   'rhombicuboctahedron', 'truncated cuboctahedron', 'triangular prism',
%   'hexagonal prism', 'octagonal prism' and 'decagonal prism'.
%3) 'plot': 'selecthinges', 'info', 'result', 'savedata' or 'plot'. 
%   'selecthinges' runs the program to select all possible hinges for a
%   defined template. 'info' should be used when setting up a new template,
%   while 'result' will run the simulation and fold the structure saving 
%   the results in .mat files. 'savedata' will read the results, save 
%   the analysis on .csv files and can plot the defomrmation of the 
%   structure. 'plot' only does the latter.
%
%See the intiOpt() file for a few other (less important) options. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GEOMETRY OPTIONS DEFAULT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input type that defines from which database to choose from
%Name of the template
opt.template='cube';
%Turn periodic boundary conditions on or off
%Options: 'on' or 'off'
opt.periodic='off';
%Extrusion length
%Value >0
opt.Lextrude=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOT AND MOVIE OPTIONS DEFAULT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type of results to display         
%Options: 'info' or 'result'
opt.plot='result';
%Number of frames to show the modes the modes
%Integer value > 0
opt.frames=30;
%Scales the amplitude of the modes
%Value between 0 and 1: 0 no deformation, 1 some hinge angles will be zero
opt.scale=0.75;
%Save all results in a gif file
%Options: 'on' or 'off'
opt.saveMovie='off';
%When using periodic boundary condtions, number of unit cells to show
%Integer value >0
opt.plotPer=2;
%Plot options default view
opt.AZ=-54;
opt.EL=14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ANALYSIS OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%UPDATE OPTIONS BASED ON INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:2:nargin
    opt.(varargin{i}) = varargin{i+1};
end
