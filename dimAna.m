%% Dimensional Analysis
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
% 
%Modified from:
%Stefan Thanheiser, "Extended Heat Transfer Model Software”. Zenodo, Aug. 
%05, 2025. doi: 10.5281/zenodo.16748359.
%
% 
%Part of the thesis:
%
%Thanheiser, S.
%A Contribution to the Development of an Active Particle Thermal Energy 
%Storage System
%PhD Thesis, TU Wien, Austria, 2025
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.17288663
%
%
%
%This script sets up the dimensional sets shown in the thesis' methodology
%section.
%
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - None


%% Influencing factors from Molerus and Wirth
h=struct('name','h','unit','W/m²K');                        %Heat transfer coefficient
c_g=struct('name','c_g','unit','J/kgK');                    %Fluidization gas isobaric specific heat capacity
c_p=struct('name','c_p','unit','J/kgK');                    %Particle isobaric specific heat capacity
rho_g=struct('name','rho_g','unit','kg/m³');                %Fluidization gas density
rho_p_rho_g=struct('name','rho_p_rho_g','unit','kg/m³');    %Difference between particle and fluidization gas density
my_g=struct('name','my_g','unit','Pas');                    %Fluidization gas dynamic viscosity
k_g=struct('name','k_g','unit','W/mK');                     %Fluidization gas thermal conductivity
w_e=struct('name','w_e','unit','m/s');                      %Excess fluidization velocity
w_mf=struct('name','w_mf','unit','m/s');                    %Minimum fluidization velocity
eps_mf=struct('name','eps_mf','unit','');                   %Minimum fluidization bed voidage (actually: 1-eps_mf)
g=struct('name','g','unit','m/s²');                         %Gravitational acceleration


%% Simple set
% D=I

%Dimensional matrices
[A,B,C,D]=dimMat(h,c_p,c_g,rho_g,w_e,w_mf,eps_mf,...
            k_g,rho_p_rho_g,g,my_g);


%Dimensional set
setSimple=dimSet(A,B,C,D);


%% Molerus and Wirth's original set
% Modifications compared to the simple set:
% Molerus and Wirth's pi_2=pi_2^-1
% Molerus and Wirth's pi_6=replace pi_6 with pi_5/pi_6
% Molerus and Wirth's pi_5=replace pi_5 with pi_5*pi_2^(1/3)
% or: small variations of the D-matrix (=linear combinations of pi-factors)

%Dimensional matrices
[A,B,C,D]=dimMat(h,c_p,c_g,rho_g,w_e,w_mf,eps_mf,...
            k_g,rho_p_rho_g,g,my_g);


%Modifications of D
D{'pi2','c_p'}=-1;
D{'pi5','c_p'}=1/3;
D{'pi6','w_e'}=1;
D{'pi6','w_mf'}=-1;


%Dimensional set
[setMolerus,C]=dimSet(A,B,C,D);


%% Auxiliary functions
%This function creates the matrices that make up the dimensional set
function [A,B,C,D]=dimMat(varargin)
    %Auxiliary functions to count dimensions from the created unit string
    count=@(str,dim) numel(regexp(str,['*?[^/]',dim]))-...
                    numel(regexp(str,['/',dim]));

    dimsum=@(nomden,dim) count(nomden{1},dim)-count(nomden{2},dim);


    %Get variables
    vars=[varargin{:}];
    n=numel(vars);

    
    %Set up dimensional matrix
    dimnames={'m','kg','s','K'};
    mat=table('Size',[length(dimnames),n],...
                'VariableTypes',repmat({'double'},1,n),...
                'VariableNames',{vars.name},...
                'RowNames',dimnames);


    %Extract dimensions from every unit
    for i=1:n
        %Normalize exponents
        vars(i).unit=strrep(vars(i).unit,'²','^2');
        vars(i).unit=strrep(vars(i).unit,'³','^3');


        %Replace derived SI units with basic SI units
        vars(i).unit=strrep(vars(i).unit,'Pa','N*m^-2');
        vars(i).unit=strrep(vars(i).unit,'W','J*s^-1');

        vars(i).unit=strrep(vars(i).unit,'J','N*m');

        vars(i).unit=strrep(vars(i).unit,'N','kg*m*s^-2');


        %Split into nominator and denominator
        nomden=strsplit(vars(i).unit,'/');
        if isscalar(nomden)
            nomden=[nomden,{''}]; %#ok<AGROW>
        end


        %Replace exponents with string of dimensions, normalize result
        for j=1:length(nomden)
            %Positive exponents
            nomden{j}=regexprep(nomden{j},...
                        '(\w+)\^(\d+)',...
                        '${repmat([$1,''*''],1,str2num($2))}');

            %Negative exponents
            nomden{j}=regexprep(nomden{j},...
                        '(\w+)\^-(\d+)',...
                        '${repmat([''/'',$1],1,str2num($2))}');

            %Cleanup and normalize
            nomden{j}=strrep(nomden{j},'*/','/');
            nomden{j}=['*',nomden{j}];
        end


        %Count units
        for j=1:length(dimnames)
            mat{dimnames{j},i}=dimsum(nomden,dimnames{j});
        end
    end
    

    %Remove empty dimensions
    mat(all(mat{:,:}==0,2),:)=[];


    %Constituting numbers
    nDims=height(mat);  %Number of dimensions
    nB=n-nDims;         %Number of variables in Matrix B
    % nA=nDims;         %Number of variables in Matrix A
    % nP=nB;            %Number of pi-factors
    
    
    %Matrices
    A=mat(:,nB+1:end);
    B=mat(:,1:nB);

    pinames=compose('pi%d',1:nB);
    C=table('Size',[nB,nDims],...
            'VariableTypes',repmat({'double'},1,nDims),...
            'VariableNames',A.Properties.VariableNames,...
            'RowNames',pinames);
    
    D=table('Size',[nB,nB],...
            'VariableTypes',repmat({'double'},1,nB),...
            'VariableNames',B.Properties.VariableNames,...
            'RowNames',pinames);
    
    D{:,:}=eye(nB);
end


%This function calculates the C matrix and creates the dimensional set
function [set,C]=dimSet(A,B,C,D)
    %Fundamental equation
    Cmat=-D{:,:}*(A{:,:}^-1*B{:,:})';
    
    
    %Fix rounding issues
    Cmat(abs(Cmat)<1e-6)=0;
    C{:,:}=Cmat;
    
    
    %Dimensional set
    set=[B,A;D,C];
end




