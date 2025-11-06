
function [fields,snaptimes,vels,CL ] = loadSnapshotFVM(filename)
    %% load FVM snapshot data
    % also returns a "cell object" structure that can be used for
    % processing
    % later as if it was experimental data

  % data = dlmread(filename);
   %% updated to use readtable instead of dlmread  
   data = readtable(filename,'ReadVariableNames', false,'Delimiter',' ','MultipleDelimsAsOne', true,'LeadingDelimitersRule', 'ignore');
   %%
   nfield = data{1,1};
   ncell = data{1,2};
   maxdeg = data{1,3};
   
   % are velocities included in these snapshots as well?
   readvel = data{1,5};
   
   % number of lines in each snapshot
   if (readvel>0) 
        snaplines = nfield+2;
   else
       snaplines = nfield+1;
   end
    
   nsnap = size(data,1)/snaplines;
  %% 
   snaptimes = data{1:snaplines:end,4};   
   if (readvel>0)
        vels = data{2:snaplines:end,:};
   else
       vels = NaN;
   end
   
   %%
   fields = zeros(nsnap,ncell,nfield);
   if (readvel>0) 
       start = 2;
   else
       start = 1;
   end
   for fc = 1:nfield
        fields(:,:,fc) = data{start+fc:snaplines:end,1:ncell};
   end
   
   % reorder dimensions as: cell, field, time
   fields = permute(fields,[2,3,1]);
   
   if (readvel>0)
       % reorganize velocities into dimensions: cell, boundary, time
       vels = reshape(vels,[nsnap,ncell,maxdeg]);
       vels = permute(vels,[2,3,1]);
   end
   
   % create a "Cell object" structure containing some minimal fields
   if (nargout>3)
       CL = CellObjPA(filename);       
       CL.NFrame = nsnap;       
       CL.dt = snaptimes(2)-snaptimes(1);
       CL.startPA = 1;
   end
end