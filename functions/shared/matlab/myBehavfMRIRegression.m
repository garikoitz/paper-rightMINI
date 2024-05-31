function myBehavfMRIRegression(trt, subs, subject_index, LD, kkvertex, tempmgh)

    % trt = 'TEST'
    % subs = TESTsubs
    % subject_index = TESTind
    
    % trt = 'RETEST'
    % subs = DAY2subs
    % subject_index = DAY2ind


% Vertexwise glmfit for PER and LEX averaged signals
labeldir      = fullfile(MINIPath, 'DATA', 'fslabeldir');

ldnames = {'CSzRT','WHzRT'};
contrastes   = {'RWvsCB','RWvsPS','RWvsSD','RWvsCS','RWvsFF','RWvsPW','RWvsNull'};
designs = {'block'};
fMRIareas = {'VOT'};
orig_state = warning('off','all');
fshome = '/Applications/freesurfer';
fssubdir = fullfile(fshome, 'subjects');

setenv('FREESURFER_HOME', fshome);
setenv('SUBJECTS_DIR', fssubdir);

brain = {};
persembrain = {};
tipocon = {'PERCEPTUAL','SEMANTIC'};  % 'PERCEPTUAL','SEMANTIC','RWvsNull'

% Create the intermediate data
for area = fMRIareas; for design = designs; for contr=1:length(contrastes)
    todos    = [subs.([ area{:} '_' design{:} '_' contrastes{contr}])]';
    todos    = todos(:, kkvertex.(area{:}));
    brain{contr}    = repmat(tempmgh.vol, [size(todos,1),1]);
    brain{contr}(:, kkvertex.(area{:})) = todos;
    if contr == 3
        % Lo tengo que hacer por sujeto independientemente
        persembrain{1} = repmat(tempmgh.vol, [size(todos,1),1]);
        for ii=1:size(brain{contr},1)
            persembrain{1}(ii,:) = sum([brain{1}(ii,:);brain{2}(ii,:);brain{3}(ii,:)], ...
                                      1,'omitnan');
        end
    end
    if contr == 6
        % Lo tengo que hacer por sujeto independientemente
        persembrain{2} = repmat(tempmgh.vol, [size(todos,1),1]);
        for ii=1:size(brain{contr},1)
            persembrain{2}(ii,:) = sum([brain{4}(ii,:);brain{5}(ii,:);brain{6}(ii,:)], ...
                                       1,'omitnan');
        end
    end
    
    % if contr == 7; persembrain{1} = brain{7}; end;        
end;end;end;




% Now that we have the averages, do the glmfit at the vertex level
thp = '13';  % '13', '20'
thpcomma = '1.3';  % '1.3', '2.0'
cwpvalthresh = '0.05';  % '0.05', '0.01'
sig = 'abs';  % 'pos', 'abs', 'neg'
ver = 'v01';
% Hacer el calculo glm
for ldns = 1:length(ldnames)
    behavname = ldnames{ldns};
    y = LD.(behavname);
    y = y(subject_index);
    ResultFldr = fullfile(MINIPath, 'DATA', 'LDfMRI', ['glmfit_' ver '_VOT_block_zSum_vertexclusterCor' thp '_' behavname]);
    if ~exist([ResultFldr]);mkdir([ResultFldr]),end
    for psb = 1:size(persembrain,2)
        tipoclust = tipocon{psb};
        Fp     = zeros(size(tempmgh.vol));
        for kk = kkvertex.(area{:})'  % Matlab = kk, FREEVIEW = kk-1
            temp1 = fitlm([persembrain{psb}(:,kk)], ... %,x3(:,kk)], ...
                          y, ...
                          'linear'); %, ... % 'linear', 'interactions'
            p_signOfRelationship = 1;
            if temp1.Coefficients.Estimate(2) > 0 
                p_signOfRelationship = -1;
            end
            Fp(kk) = p_signOfRelationship * log10(temp1.coefTest);
        end
       % Write mgh 
       logP = tempmgh;
       logP.vol = Fp;
       inFile = [trt '_' behavname  '_' tipoclust];
       MRIwrite(logP, [ResultFldr fsp inFile '.mgh']);
end;end
        
% Montecarlo (do it once and store it for later use)
%   First create the simulation inside our ROI. 
% It was run in the server and then copied locally, created with the following
% command line
% correc = ['mri_mcsim --o $FREESURFER_HOME/average/mult-comp-cor/fsaverage/lh/VOT ' ...
%           '--base mc-z --surface fsaverage lh --nreps 10000 ' ...
%           '--label  VotNoV1V2yMin16'];
      



for ldns = 1:length(ldnames)
    behavname = ldnames{ldns};
    for psb = 1:size(persembrain,2)
        tipoclust = tipocon{psb};
        inFile = [trt '_' behavname  '_' tipoclust];
        ResultFldr = fullfile(MINIPath, 'DATA', 'LDfMRI', ['glmfit_' ver '_VOT_block_zSum_vertexclusterCor' thp '_' behavname]);
        cd(ResultFldr); if ~exist(inFile); mkdir(inFile); end

        cmdsc = [fsbin fsp 'mri_surfcluster ' ...
             '--in ' inFile '.mgh ' ...
             '--csd ' fullfile(MINIPath,'DATA','MCsims','lh',area{:},'fwhm05',sig,['th' thp],'mc-z.csd ') ...
             '--mask ' fullfile(MINIPath,'DATA','MCsims','lh', area{:}, 'mask.mgh ') ...
             '--cwsig ' inFile fsp 'cache.th' thp '.' sig '.sig.cluster.mgh ' ...
             '--vwsig ' inFile fsp 'cache.th' thp '.' sig '.sig.voxel.mgh ' ...
             '--sum ' inFile fsp 'cache.th' thp '.' sig '.sig.cluster.summary ' ...
             '--ocn ' inFile fsp 'cache.th' thp '.' sig '.sig.ocn.mgh ' ...
             '--oannot ' inFile fsp 'cache.th' thp '.' sig '.sig.ocn.annot ' ...
             '--annot aparc ' ...
             '--csdpdf ' inFile fsp 'cache.th' thp '.' sig '.pdf.dat ' ...
             '--cwpvalthresh ' cwpvalthresh ' ' ...
             '--thmin ' thpcomma ' ' ...
             '--o ' inFile fsp 'cache.th' thp '.' sig '.sig.masked.mgh ' ...
             '--no-fixmni ' ...
             '--bonferroni 2 ' ...
             '--surf white ']
         system(cmdsc)

    end
end
% Extract the labels for the annotation
for ldns = 1:length(ldnames)
    behavname = ldnames{ldns};
    for psb = 1:size(persembrain,2)
        tipoclust = tipocon{psb};
        inFile = [trt '_' behavname  '_' tipoclust];
        ResultFldr = fullfile(MINIPath, 'DATA', 'LDfMRI', ['glmfit_' ver '_VOT_block_zSum_vertexclusterCor' thp '_' behavname]);
        cd(ResultFldr); cd(inFile); 
        cmdsc = [fsbin fsp 'mri_annotation2label ' ...
                    '--subject fsaverage ' ...
                    '--hemi lh ' ...
                    '--annotation ' fullfile(ResultFldr,inFile,['cache.th' thp '.' sig '.sig.ocn.annot ']) ...
                    '--outdir ' fullfile(ResultFldr, inFile)];

        system(cmdsc)
    end
end


%% Ahora hacer las mismas correlaciones 2 (SEM,PER) x 2 (CS, RW) dentro de los ROIs perVWFA y semVWFA, y dentro de los clusters principales
thp = '13';  % '13', '20'
thpcomma = '1.3';  % '1.3', '2.0'
cwpvalthresh = '0.05';  % '0.05', '0.01'
sig = 'abs';  % 'pos', 'abs', 'neg'
ver = 'v05';
sdOutlier = 2;
tipocon = {'PERCEPTUAL','LEXICAL'};
% regions = {'semVWFA4', 'perVWFA4'};  
regions = {[trt '_PER_CS'], [trt '_PER_RW'],[trt '_LEX_CS'], [trt '_LEX_RW']};  
ldnames = {'CSzRT','WHzRT'};
% myFig = figure();



TableElements = {'SUBJECT','LDcat','fMRIcat','LD','spmT'};
datos = array2table(zeros(0,length(TableElements)));
datos.Properties.VariableNames = TableElements;

ii = 0;
for ldns = 1:length(ldnames)
    behavname = ldnames{ldns};
    y = LD.(behavname);
    y = y(subject_index);
    % Clean outliers in behav data
    % sum(isnan(y));
    % lowlim  = (mean(y,'omitnan') - sdOutlier*std(y,'omitnan'));
    % highlim = (mean(y,'omitnan') + sdOutlier*std(y,'omitnan'));
    % y(y < y | y > y) = NaN;
    % sum(isnan(y));
    for psb = 1:size(persembrain,2)
        tipoclust = tipocon{psb};
        for reg = regions
            ii = ii+1;
            bbrain = persembrain{psb}(:,kkvertex.(reg{:}));
            if isequal(ones(1,1), kkvertex.(reg{:})); bbrain(bbrain==0)=NaN;end
            mbbrain    = mean(bbrain, 2, 'omitnan');
            % Sabemos que los sujetos 13 y 18 para block son NaN
            mbbrain(13) = NaN;
            mbbrain(18) = NaN;
            % ver por si acaso si hay outliers o no:
            % en behav
            sum(isnan(mbbrain));
            lowlim  = (mean(mbbrain,'omitnan') - sdOutlier*std(mbbrain,'omitnan'));
            highlim = (mean(mbbrain,'omitnan') + sdOutlier*std(mbbrain,'omitnan'));
            mbbrain(mbbrain < mbbrain | mbbrain > mbbrain) = NaN;
            sum(isnan(mbbrain));
            [resultado,pval] = corr(mbbrain, y,'rows','pairwise');
            disp([behavname '_' tipoclust '_' reg{:} ': corr=' num2str(resultado) ', pValue=' num2str(pval)]);
            % min(mbbrain)
            % max(mbbrain)
            
            
            listToPlot = [7,16,1,10];
            if any(ii==listToPlot)
%                 donde =   [2  ,0  ,0  ,0  ,0  ,0  ,1  ,0  ,0  ,4  ,0  ,0  ,0  ,0  ,0  ,3];
%                 koloria = ['b','g','g','g','g','g','r','g','g','b','g','g','g','g','g','r'];
%                 xticks = {};
%                 xticks{7} = [-2:2:4]; xticks{16} = [-2:2:4];
%                 xticks{1} = [-4:2:10]; xticks{10} = [-4:2:10];
%                 xlims = {};
%                 xlims{7} = [-3,4.5]; xlims{16} = [-3,4.5];
%                 xlims{1} = [-4.3,10]; xlims{10} = [-4.3,10];
                
                tmpsubs = struct2cell(subs);
                tmpsubs = tmpsubs(1,:)';
                tmpdatos = cell2table(tmpsubs);
                tmpdatos.Properties.VariableNames = {'SUBJECT'};
                LDcat = behavname(1:2); if strcmp(LDcat,'WH');LDcat = 'RW';end
                tmpdatos.LDcat   = repmat(LDcat, [height(tmpdatos),1]);
                tmpdatos.fMRIcat = repmat(tipoclust(1:3), [height(tmpdatos),1]);
                tmpdatos.LD      = y;
                tmpdatos.spmT    = mbbrain;
                % Add it to the table
                datos  = [datos; tmpdatos];
                
%                subplot(2,2,donde(ii));
%                     s1 = plot(mbbrain,y,'o','color',[0 0 0],'markerfacecolor',koloria(ii),'markersize',6);
%                     h = lsline; set(h(1),'color', koloria(ii), 'LineWidth', 2, 'MarkerSize', 15);
%                     title({['ROI ' strrep(reg{:},'_','\_')], ...
%                            [behavname '_' tipoclust '_' reg{:} ': corr=' num2str(resultado) ', pValue=' num2str(pval)]});
%                     xlabel(['fMRI ' tipoclust]);
%                     ylabel(['RT ' behavname]);
%                     % axis scaling
%                     set(gca, 'YTick',[-2:1:1]);
%                     ylim([-2.5 1.5]); 
%                     set(gca,'XTick', xticks{ii});
%                     xlim(xlims{ii});
            end
        end
    end
end
% Write the files to be used in R, plot there...
writetable(datos, ...
    fullfile(MINIPath,'DATA','LDfMRI',[trt '_Data4BehavfMRIRegression.csv']), ...
    'FileType', 'text', ...
    'Delimiter', 'comma', ...
    'WriteVariableNames', true); 

    
    
    
    
    

end