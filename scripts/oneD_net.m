%line.net
addpath("~/ER_Ca_rot/networkSpreadFVM-main/scripts/");
addpath("~/ER_Ca_rot/networkSpreadFVM-main/examples/");


nodelblfmtstring = ['NODE %d ' repmat(['%20.10f '],1,3) '%s \n'];
nodefmtstring = ['NODE %d ' repmat(['%20.10f '],1,3) '\n'];
edgefmtstring = ['EDGE %d %d %d %20.10f\n'];
for k = 10:10:150 
    net_file = "~/ER_Ca_rot/networkFVMsims_git/examples/param_1D/line"+k+".net";
    of = fopen(net_file,'w');
     nnodes = k;
    pos=zeros(nnodes,3);%xpos
    edges=zeros(nnodes-1,3);
   
    for i = 1:k
         pos(i,1) = (i-1);
        pos(i,2) = 0.5;
        if(i==18)
            fprintf(of,nodelblfmtstring,(i), pos(i,:),'FN');
        else
            fprintf(of,nodefmtstring,(i), pos(i,:));
    
        end
    
    end
for n = 1:k-1
    
    fprintf(of,edgefmtstring,(n), [n,n+1,1.00]);
    
end

end
fclose(of);
% NT.nodepos = pos/40;
% NT.nnode = nnodes;
% nodel = strings(NT.nnode,1);
% nodel(1) ="FN";
% NT.nodelabels = nodel;
% NT.setupNetwork();
% NT
% NT.outputNetwork(net_file);
%d = NT.dim;

%%
10^-6*6*10^23*40*10^-6*3*2500*10^-18
