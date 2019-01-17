%modified Marsbar Image scaling (SPM12 compatible)
%imported from Marseille connectivity project
% 17-02-2017 Jonathan Wirsich / Connectlab
function [VY,row] = SPMmars_image_scaling(marsD, subjectpath, filepattern)
% get image scaling data for images, maybe via SPM design
% FORMAT [VY,row] = mars_image_scaling(marsD)
%-----------------------------------------------------------------------
%
% Inputs
% marsD      - design matrix to (optionally) get parameters from
% 
% Returns
% VY         - SPM vol structs with selected scaling
% row        - cell array, one per subject/session giving corresponding
%              rows in for VY array
%
% $Id$

VY = [];
if nargin < 1
  marsD = [];
end

%fMRI
dGM =   100;
sess_str = 'Session';
  
VY = [];
Global = [];  

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Extract data from ROI(s)');

% images
% if isempty(marsD)
  spmf = 0;
% else
%   spmf = spm_input('Images from:', '+1','b',['SPM design|GUI select'], ...
% 		   [1 0], 2);
% end

% get images, from design, or by hand
% if spmf
%   if ~has_images(marsD);
%     warning('Design structure does not specify images');
%     return
%   end
%   VY = get_images(marsD);
%   row = block_rows(marsD);
%   nsess = length(row);
% end

if isempty(VY)  % need to know about images
  % no of sessions / subjects
  nsess = 1;
%   nsess = spm_input(sprintf('No of %ss', sess_str), '+1', 'r', 1, 1); 
  % select files for each session
  for s = 1:nsess
    simgs = spm_select('FPlist', subjectpath, filepattern);
    row{s} = (1:size(simgs, 1))'+size(VY,1);
    VY = strvcat(VY, simgs);
  end 
end  % of image get routines
if isempty(VY), return, end

% global scaling options
askGMf = 1;
% scaling by hand
%   glabs = [sess_str ' specific scaling',...
% 	   '|Proportional scaling|Raw data'];
%   tmp = spm_input('Scaling from:', '+1', 'm', glabs, [1 2 3], 1);

    Global = [];

% Grand mean scaling
GM = 0;
% if askGMf
%   GM = spm_input('Scale grand mean to (0=raw)','+1','r',dGM,1);
% end

% map files now, if not yet mapped
if ischar(VY)
  fprintf('\n%-40s: %30s','Mapping files',' ')                     %-#
  VY = spm_vol(VY);
  fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')         %-#
end

% Apply scaling options if necessary
if ~isempty(Global)  
  
%-Compute Global variate
%-----------------------------------------------------------------------
q      = length(VY);
g      = zeros(q,1);
fprintf('%-40s: %30s','Calculating globals',' ')                     %-#
for i  = 1:q
  fprintf('%s%30s',repmat(sprintf('\b'),1,30),sprintf('%4d/%-4d',i,q)) %-#
  g(i) = spm_global(VY(i));
end
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')               %-#

% get null GM scaling
if (GM == 0)
  GM = mean(g);
end

% scale if specified (otherwise subject / session specific grand mean scaling)
%-----------------------------------------------------------------------
gSF     = GM./g;
if strcmp(Global,'None')
  for i = 1:nsess
    j      = row{i};
    gSF(j) = GM./mean(g(j));
  end
end

%-Apply gSF to memory-mapped scalefactors to implement scaling
%-----------------------------------------------------------------------
for  i = 1:q, VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*gSF(i); end

end % of global options


  