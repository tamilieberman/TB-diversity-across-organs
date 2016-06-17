close all;


%% Parameter set up


%for genotype identification
min_snp_freq_in_single_sample_for_genotype_placement=.2;
purity_cutoff_for_sample_sorting=.8;
purity_cutoff_for_adding_in_singletons = .8 ;
min_average_coverage_for_genotype_id = 12;
loose_freq_threshold_for_majority = .4;
strict_freq_threshold_for_majority = .6; %strict treshold to avoid false groupings
noise_buffer_for_singletons =.15 ;

min_freq_genotype_in_at_least_one_site = .02;

%for genotype assignment
snp_threshold_for_assigning_to_genotypes = .2;
max_remaining_room_to_add_in_ancestor = .15;
max_unexplained_snps_to_add_in_ancestor = .8;

%for mixed strain infection plotting and analysis
number_dummy_muts_for_mixed_infections = 20;

%for transmission analysis
min_observations_for_transmission_analysis = 2;
min_isolates_per_site = 7;
min_freq_within_site = .03;

% how far upstream of the nearest gene to annotate something a promoter mutation
promotersize=150;


%% Enviornment set up -- probably won't need to change

masterdir=char(pwd);

REFGENOMEFOLDER=[masterdir '/MTB_anc'];
SCRIPTSDIRECTORY = [masterdir '/scripts'];
path(SCRIPTSDIRECTORY,path);


%% Things I need to convert into csv files

subjects=1:44;
subjects_with_multiple_strains=[11 15 28 37];

sites={'Lung1','Lung2','Lung3','Lung4','Lung5','Lung6','Eta','Liver','Spleen','SerousFluid','SerousFluid2','Lymph','LungPus'};
site_colors={'MidnightBlue','CornflowerBlue', 'Blue', 'DarkTurquoise','LightBlue','Teal','Green','DeepPink','Purple','Orange','DeepPink', 'Yellow', 'Maroon','Maroon', 'Maroon', 'Maroon', 'Maroon','Maroon'};


%% Initialize

NTs='ATCG';


folders={'trees' 'trees/genotypetrees_lung', 'trees/genotypetrees_organs', 'trees/temptreefolder', 'trees/all', 'summary_graphics', 'csvfiles_generated'};
for i=1:numel(folders)
    if ~exist(folders{i},'dir')
        mkdir(folders{i})
    end
end

if ~exist('trees/dnapars','file')
    copyfile([SCRIPTSDIRECTORY '/dnapars'],'trees/dnapars')
end

eta_site_number=find(strcmp(sites,'Eta'));
liver_site_number=find(strcmp(sites,'Liver'));
spleen_site_number=find(strcmp(sites,'Spleen'));


fid_tableS1=fopen([masterdir '/csvfiles_generated/forsupptable2.csv'],'w');
fid_tableS2=fopen([masterdir '/csvfiles_generated/forsupptable3.csv'],'w');
fid_tableS3=fopen([masterdir '/csvfiles_generated/forsupptable4.csv'],'w');


%% Make empty data structures

observed_matrices={};
de_novo_muts_all=[];
de_novo_muts_placed=[];
adjacent_removed=[];

genotypematrices={};
genotypes_freq_sites_matrices={};
genotype_sites_times_matrices={};
compartment_lists={};

number_mutations_seperating_strains={};
distance_to_subject_MRCA=zeros(length(subjects),1);
observed_matrices_placed={};
observed_matrices_placed_with_multiple_lineages={};
average_coverage_isolates={};
defined_lineages_by_subject={};
genotypematrices_all_mutations={};
proportionlineage1_by_subject={};


genes_mutated_multiple_times_within_a_subject=[];
operons_mutated_multiple_times_within_a_subject=[];

subject_lineages=zeros(numel(subjects),1);
samples_per_subject=zeros(numel(subjects),1);


proprotion_minor_lineage=zeros(length(subjects),1);
average_coverage_subject=zeros(length(subjects),1);
number_specimens=zeros(length(subjects),1);

num_duplicates_removed=zeros(numel(subjects),1);
num_removed_because_name_not_well_formed=zeros(numel(subjects),1);

max_distance_between_genotypes=zeros(numel(subjects),1);
max_distance_between_genotypes_in_eta=nan(numel(subjects),1);


%% Call genotypes by subject

for k=1:numel(subjects)
    
    
    %% load in
    cd([masterdir '/subject_folders/P' num2str(k)])
    load('de_novo_muts')
    observedmatrix(isnan(observedmatrix))=0;
    numsamples=length(SampleNames);
    
    %% categorize samples by locations
    
    isolate_compartments=zeros(numel(SampleNames),1);
    for i=1:numel(SampleNames)
        splitname=strsplit(SampleNames{i},'-');
        isolate_compartments(i)=find(strcmp(splitname{2},sites));
    end
    
    %% For Supplementary table 2
    
    for i=1:spleen_site_number
        if sum(isolate_compartments==i) > 0
            fprintf(fid_tableS2,[num2str(sum(isolate_compartments==i)) ',']);
        else
            fprintf(fid_tableS2,',');
        end
    end
    
    for i=(spleen_site_number+1):numel(sites)
        if sum(isolate_compartments==i) > 0
            fprintf(fid_tableS2,[num2str(sum(isolate_compartments==i)) ',' sites{i} ',']);
        end
    end
    fprintf(fid_tableS2,'\n');
    
    
    
    
    %% save basic stats
    
    average_coverage_subject(k)=mean(averagecoverage);
    samples_per_subject(k)=numel(SampleNames);
    number_specimens(k)=numel(unique(isolate_compartments));
    
    if ~isempty(de_novo_muts)
        
        
        %Sanity check for adjacent nucleotides
        gp=[de_novo_muts(:).pos];
        adjacent_starts=find(gp(2:end)-gp(1:end-1) < 2);
        if ~isempty(adjacent_starts)
            error('Error: Directly adjacent mutations found in single patient. Probably includes false positives. ')
        end
        
        %add this field to make things easier later
        for i=1:numel(de_novo_muts)
            de_novo_muts(i).subject=k;
        end
        
        %% Add in lineage info
        if mean(proportionlineage1)~=1
            observedmatrix=[observedmatrix; repmat(proportionlineage1,20,1); repmat(1-proportionlineage1,20,1)];
        end
        
        
        %% Find genotypes
        
        [genotypematrix, sample_that_defined_each_genotype,orderbycoverage] = find_genotypes(observedmatrix, averagecoverage, ismember(k,subjects_with_multiple_strains), number_dummy_muts_for_mixed_infections,...
            min_snp_freq_in_single_sample_for_genotype_placement, purity_cutoff_for_sample_sorting, purity_cutoff_for_adding_in_singletons,...
            min_average_coverage_for_genotype_id,loose_freq_threshold_for_majority, strict_freq_threshold_for_majority, noise_buffer_for_singletons);
        
        
        %% Assign genotypes to samples
        
        [genotypes_freq_sites, genotype_sites_times, calls_for_tree, TreeSampleNames, genotypeOrder] = ...
            assign_genotypes(observedmatrix,isolate_compartments,...
            genotypematrix,numel(sites),snp_threshold_for_assigning_to_genotypes,...
            max_remaining_room_to_add_in_ancestor,max_unexplained_snps_to_add_in_ancestor, SampleNames);
        
        
        %% Make parsimony tree with all called genotypes called in all samples -- Figure 2C
        
        nts_for_tree=char(calls_for_tree);
        nts_for_tree(calls_for_tree==0)='A';  %actually identity doesn't matter -- 0's and 1's is fine
        nts_for_tree(calls_for_tree==1)='T';
        
        cd([masterdir '/trees/temptreefolder'])
        treefilename=generate_parsimony_tree(nts_for_tree, TreeSampleNames, ['P' num2str(k)]);
        TreeColors={}; for i=1:numel(TreeSampleNames); TreeColors{i}='16777216';end;
        color_tree(treefilename, [masterdir '/trees/all/P' num2str(k) '.tree'], TreeSampleNames, TreeColors)
        cd([masterdir '/subject_folders/P' num2str(k)])
        
        
        
        
        %% Visualize each subject's mutation matrix
        
        h=figure(6); clf; hold on;
        
        %Sort samples and show their names
        sites_in_subject=unique(isolate_compartments);
        [sortedsites, orderisolates]=sort(isolate_compartments);
        for i=1:numel(sites_in_subject)
            c=sites_in_subject(i);
            plot([-3 numel(de_novo_muts)+1], [find(sortedsites==c,1) find(sortedsites==c,1)]-.5 ,'k-')
            numsamplesc=sum(isolate_compartments==c);
            text(-3,find(sortedsites==c,1)+(numsamplesc/2 -.5),sites{c}, 'FontSize', 12)
        end
        plot([-3 numel(de_novo_muts)+1], [numel(SampleNames) numel(SampleNames)]+.5 ,'k-')
        
        
        %Plot information about multiple strain infection and global
        %lineages
        if mean(proportionlineage1)==1 %single lineages
            title(['P' num2str(k) '   - Lineage ' num2str(find(mean(proportion_global_lineages(1:6,:),2)>.95))],'FontSize',12)
            subject_lineages(k)=find(mean(proportion_global_lineages(1:6,:),2)>.95);
        else
            lineages=find(mean(proportion_global_lineages(1:6,:),2)>.05);
            t=['P' num2str(k)  '   - mixed Lineage ' num2str(lineages(1))];
            if numel(lineages)>1
                t=[t ' &  Lineage ' num2str(lineages(2))];
            end
            title(t,'FontSize',12)
            for s=1:numel(proportionlineage1)
                y=find(orderisolates==s);
                proportionlineage1(isnan(proportionlineage1))=0;
                patch([-2 -2 -1 -1],[y-.5 y+.5 y+.5 y-.5],'r', 'FaceAlpha', proportionlineage1(s), 'FaceColor', rgb('Red'))
                patch([-1 -1 0 0],[y-.5 y+.5 y+.5 y-.5],'r', 'FaceAlpha', 1-proportionlineage1(s), 'FaceColor', rgb('Blue'))
            end
        end
        
        for i=1:numel(de_novo_muts)
            %Plot information about the location of each mutation
            if isempty(de_novo_muts(i).locustag)
                positionlabel=[de_novo_muts(i).gene1 '-' de_novo_muts(i).gene2]; %Intragenic
            else
                if ~isempty(de_novo_muts(i).gene)
                    genename=de_novo_muts(i).gene;
                else
                    genename=de_novo_muts(i).locustag;
                end
                if ~isempty(de_novo_muts(i).muts)  %Nonsynonymous
                    positionlabel=[genename '-' de_novo_muts(i).muts{1}];
                else %Synonymous
                    positionlabel=[genename '(S)'];
                end
            end
            
            %Show mutation frequency matrix
            samples_p=find(observedmatrix(i,:)>.08);
            for s=samples_p
                y=find(orderisolates==s);
                patch([i i i+1 i+1],[y-.5 y+.5 y+.5 y-.5],'r', 'FaceAlpha', observedmatrix(i,s), 'FaceColor', rgb('Black'))
            end
            text(i+.5,-.15*(numel(SampleNames)),positionlabel, 'FontSize', 6, 'Color',rgb('Black'), 'Rotation', 90)
        end
        
        %Plot information about which samples defined which genotypes
        for i=1:numel(sample_that_defined_each_genotype)
            if sample_that_defined_each_genotype(i) > 0  %important for multiple strain infections
                y=find(orderisolates==sample_that_defined_each_genotype(i));
                snps_i=find(genotypematrix(i,:)>0);
                for j=snps_i
                    patch([j j j+1 j+1],[y-.5 y+.5 y+.5 y-.5],'r', 'FaceAlpha', 0, 'FaceColor', rgb('Black'), 'EdgeColor', rgb('Red'), 'LineWidth', 2)
                end
            end
        end
        
        set(gca, 'Ytick',[]); xlim([-3 numel(de_novo_muts)+1]); ylim([-.15*(numel(SampleNames)) numel(SampleNames)+1])
        saveas(h,[masterdir '/summary_graphics/P' num2str(k) 'summary.jpg'],'jpg');close(h)
        
        
        
        %% make a supplementary figure 1 illustrating genotype calling in subject 20
        
        if k==20;
            
            colors_for_fig2={'Thistle','FireBrick','DeepPink','Red','Coral','Wheat','Goldenrod','Gold'...
                'MediumVioletRed','DeepPink','FireBrick','Red','Coral','Orange','SandyBrown','Yellow','Goldenrod'...
                'MediumVioletRed','DeepPink','FireBrick','Red','Coral','Orange','SandyBrown','Yellow','Goldenrod'};
            
            figure; clf; hold on;
            title(['P' num2str(k)])
            
            j=0;
            for s=1:numel(orderbycoverage)
                y=numel(orderbycoverage)-s+.5;
                if ismember(orderbycoverage(s),sample_that_defined_each_genotype)
                    defineslineage=1;
                    [~,lineagedefined]=ismember(orderbycoverage(s),sample_that_defined_each_genotype);
                    j=j+1;
                else
                    defineslineage=0;
                end
                plot([0 size(observedmatrix,1)+2],[y y],':', 'Color', rgb(site_colors{isolate_compartments(orderbycoverage(s))}),'LineWidth',2)
                for p=1:size(observedmatrix,1);
                    if observedmatrix(p,orderbycoverage(s)) >= .05 %min(.08,max(omatrix(p,:)))
                        x=p; %(1- min(omatrix(p,orderbycoverage(s))+.2,1))*[1 1 1]
                        patch([x x x+.9 x+.9],[y-.4 y+.4 y+.4 y-.4],'r','FaceColor', (1- observedmatrix(p,orderbycoverage(s)))*[1 1 1],'EdgeColor', 'none')
                        if defineslineage & genotypematrix(lineagedefined,p)>0
                            plot([x x],[y-.4 y+.4],'Color', rgb(colors_for_fig2{j}),'LineWidth',2)
                            plot([x+.9 x+.9],[y-.5 y+.4],'Color', rgb(colors_for_fig2{j}),'LineWidth',2)
                            plot([x x+.9],[y+.4 y+.4],'Color', rgb(colors_for_fig2{j}),'LineWidth',2)
                            plot([x x+.9],[y-.4 y-.4],'Color', rgb(colors_for_fig2{j}),'LineWidth',2)
                        else
                            patch([x x x+.9 x+.9],[y-.4 y+.4 y+.4 y-.4],'r', 'FaceColor', (1-observedmatrix(p,orderbycoverage(s)))*[1 1 1],'EdgeColor', 'none')
                        end
                    end
                end
            end
            set(gca,'Xtick',[])
            set(gca,'Ytick',[])
            xlim([-1 inf])
            ylim([0 numel(orderbycoverage)+1])
        end
        
        
        %%  average distance to MRCA
        distance_to_MRCA_by_site=zeros(size(sites_in_subject));
        for y=1:numel(sites_in_subject)
            distance_to_MRCA_by_site(y)=mean(sum(observedmatrix(:,isolate_compartments==sites_in_subject(y))));
        end
        
        %% max distance between genotypes
        
        maxdist = find_max_distance_between_genotypes(genotypematrix);
        if ismember(k,subjects_with_multiple_strains)
            maxdist=maxdist-(2*number_dummy_muts_for_mixed_infections)+sum(num_muts_seperating_strains);
        end
        max_distance_between_genotypes(k)=maxdist;
        
        %repeat analysis within endotracheal aspirate only
        found_in_eta=genotypes_freq_sites(:,strcmp(sites,'Eta'))>0;
        max_distance_between_genotypes_in_eta(k)=find_max_distance_between_genotypes(genotypematrix(found_in_eta(1:end-1),:));
        %don't include ancestor in case every sample in eta shares a mutation
        
        %% save everything
        compartment_lists{k}=isolate_compartments;
        observed_matrices{k}=observedmatrix;
        average_coverage_isolates{k}=averagecoverage;
        observed_matrices_only_de_novo_mutations{k}=observedmatrix;
        
        genotypematrix=genotypematrix(genotypeOrder,:); %first sort to match, because times and freqs already sorted
        genotypes_freq_sites_normalized=genotypes_freq_sites./repmat(sum(genotypes_freq_sites),size(genotypes_freq_sites,1),1); %normalize this
        
        
        %remove genotypes at very low frequency, unclassified mutations
        genotypestoremove=sum(genotypes_freq_sites_normalized(1:end-1,:)>min_freq_genotype_in_at_least_one_site,2)==0;
        mutationstoremove=sum(genotypematrix(~genotypestoremove,:),1)==0;
        
        genotypes_freq_sites(genotypestoremove,:)=[];
        genotypes_freq_sites_matrices{k}=genotypes_freq_sites;
        
        genotype_sites_times(genotypestoremove,:)=[];
        genotype_sites_times_matrices{k}=genotype_sites_times;
        
        
        
        genotypematrix(genotypestoremove,:)=[];
        genotypematrices_all_mutations{k}=genotypematrix;
        
        genotypematrices{k}=genotypematrix(:,~mutationstoremove);
        observed_matrices_placed_with_multiple_lineages{k}=observedmatrix(~mutationstoremove,:);
        
        if ismember(k,subjects_with_multiple_strains)
            mutationstoremove=mutationstoremove(1:end-2*number_dummy_muts_for_mixed_infections);
            observed_matrices_only_de_novo_mutations{k}=observedmatrix(1:end-(2*number_dummy_muts_for_mixed_infections),:);
            proportionlineage1_by_subject{end+1}=proportionlineage1;
            if mean(proportionlineage1) > .5
                proprotion_minor_lineage(k)=1-mean(proportionlineage1);
                proprotion_minor_lineage_isolates=1-proportionlineage1;
            else
                proprotion_minor_lineage(k)=mean(proportionlineage1);
                proprotion_minor_lineage_isolates=proportionlineage1;
            end
            observed_matrices_only_de_novo_mutations{k}=observedmatrix(~mutationstoremove,:);
        end
        
        
        de_novo_muts_all=[de_novo_muts_all de_novo_muts];
        de_novo_muts_placed=[de_novo_muts_placed de_novo_muts(~mutationstoremove)];
        observed_matrices_placed{k}=observedmatrix(~mutationstoremove,:);
        
        distance_to_subject_MRCA(k)=mean(distance_to_MRCA_by_site);
        if ismember(k,subjects_with_multiple_strains)
            distance_to_subject_MRCA(k)=distance_to_subject_MRCA(k)-(2*number_dummy_muts_for_mixed_infections)+mean(num_muts_seperating_strains);
        else
            defined_lineages=sample_that_defined_each_genotype(genotypeOrder);
            defined_lineages(genotypestoremove)=[];
            defined_lineages_by_subject{k}=sample_that_defined_each_genotype;
        end
        number_mutations_seperating_strains{k}=num_muts_seperating_strains;
        
        
        %% print supplementary table 1
        
        [~,l2]=max(mean(proportion_global_lineages,2));
        fprintf(fid_tableS1,[num2str(samples_per_subject(k)) ',' num2str(average_coverage_subject(k)) ',' num2str(mean(maxdist)) ','  num2str(distance_to_subject_MRCA(k)) ',' ...
            num2str(l2) ',' num2str(numel(de_novo_muts)) '\n']);
        
        
        %% print supplementary table 2
        
        fprintf(fid_tableS3,'Subject P%i\n',k);
        fprintf(fid_tableS3, 'Position,CDS,Gene Name,Annotation,Ancestral,Mutation,Coding change,');
        for i=1:numel(SampleNames)
            fprintf(fid_tableS3,[SampleNames{i} ',']);
        end
        fprintf(fid_tableS3,'\n');
        for j=1:numel(de_novo_muts)
            m=de_novo_muts(j);
            alt=m.nts; alt(alt==m.anc)=[];
            if ~isempty(m.muts)
                mut=m.muts{1};
            else
                mut='';
            end
            if ~isempty(m.locustag)
                x=m.annotation;
                x(x==',')=':';
                fprintf(fid_tableS3, '%i,%s,%s,%s,%s,%s,%s,', m.pos,m.locustag, m.gene,x, m.anc, alt, mut);
            else
                x=['Intergenic ' m.annotation];
                x(x==',')=':';
                fprintf(fid_tableS3, '%i,,,%s,%s,%s,,', m.pos,x, m.anc, alt);
            end
            for i=1:numel(SampleNames)
                fprintf(fid_tableS3,'%2.2f,', observedmatrix(j,i));
            end
            fprintf(fid_tableS3,'\n');
        end
        fprintf(fid_tableS3, ',,,,,,Coverage,');
        for i=1:numel(SampleNames)
            fprintf(fid_tableS3,'%3.1f,', averagecoverage(i));
        end
        fprintf(fid_tableS3,'\n');
        if ismember(k,subjects_with_multiple_strains)
            fprintf(fid_tableS3, ',,,,,,Proprotion minor lineage,');
            for i=1:numel(SampleNames)
                fprintf(fid_tableS3,'%2.2f,', proprotion_minor_lineage_isolates(i));
            end
        end
        fprintf(fid_tableS3,'\n\n\n');
        
        
    else
        %% print supplementary table 1
        [~,l2]=max(mean(proportion_global_lineages,2));
        fprintf(fid_tableS1,[num2str(samples_per_subject(k)) ',' num2str(average_coverage_subject(k)) ',0,0,' num2str(l2) ',0\n']);
        
        %% print supplementary table 2
        
        
        fprintf(fid_tableS3,'Subject P%i\n',k);
        fprintf(fid_tableS3, 'No de novo mutations detected,, ,,,,,');
        for i=1:numel(SampleNames)
            fprintf(fid_tableS3,[SampleNames{i} ',']);
        end
        fprintf(fid_tableS3,'\n');
        fprintf(fid_tableS3, ',,,,,,Coverage,');
        for i=1:numel(SampleNames)
            fprintf(fid_tableS3,'%3.1f,', averagecoverage(i));
        end
        fprintf(fid_tableS3,'\n');
        if ismember(k,subjects_with_multiple_strains)
            fprintf(fid_tableS3, ',,,,,,Proprotion minor lineage,');
            for i=1:numel(SampleNames)
                fprintf(fid_tableS3,'%1.2f,', proprotion_minor_lineage(i));
            end
        end
        fprintf(fid_tableS3,'\n\n\n');
        
    end
end


cd(masterdir)

%% make useful structures

genomic_positions=[de_novo_muts_placed(:).pos];
genomic_positions_all=[de_novo_muts_all(:).pos];

subject_each_pos=[de_novo_muts_all(:).subject];

%% Sanity checks

%note if the same mutation was found de novo in multiple subjects
[number_times_found,unique_genomic_positions]=find_count_duplicates(genomic_positions_all);
if sum(ismember(genomic_positions,unique_genomic_positions(number_times_found>1))) > 0
    disp('Warning: The following de novo mutations were detected in independent subjects:')
    disp(unique_genomic_positions(number_times_found>1))
end

%check that none of the inferred trees violate parsimony -- do this in a sort of exhaustive way
for k=1:numel(subjects)
    gmatrix=genotypematrices{k};
    if ismember(k,subjects_with_multiple_strains) & size(gmatrix,2)>(2*number_dummy_muts_for_mixed_infections)
        gmatrix=gmatrix(:,1:(end-(2*number_dummy_muts_for_mixed_infections)));
    end
    nonunique=find(sum(gmatrix)>1);
    if numel(nonunique)>1
        for a=1:numel(nonunique)
            %for each mutation, find each mutation it is ever in a genotype with
            i=nonunique(a);
            comuts=find(sum(gmatrix(gmatrix(:,i)>0 & sum(gmatrix,2)>0,:)) & 1:size(gmatrix,2)>i');
            for b=1:numel(comuts)
                j=comuts(b);
                if ismember([1 0],gmatrix(:,[i j]),'rows') & ismember([0 1],gmatrix(:,[i j]),'rows')
                    disp(['P' num2str(k)])
                    disp(gmatrix(sum(gmatrix(:,nonunique),2)>0,:))
                    error('Two mutations each found in seperate genotypes within the same subject (parsimony violation)')
                end
            end
        end
    end
end

%% plot distribution of strains across lung for multiple-strain infection - Fig 3A
figure; clf;
i=1;

fid_F3=fopen([masterdir '/csvfiles_generated/figure3sourcedata.csv'],'w');
fprintf(fid_F3,'#Figure 3A\n#Also see Supplementary Table 4 for Proportion Strain 1 by sample\n');
fprintf(fid_F3,'#Subject,Site,Mean(Proportion Strain 1), SEM(Proportion Strain 1)\n');
for k=1:numel(subjects)
    if ismember(k,subjects_with_multiple_strains)
        
        subplot(4,1,i);hold on;
        proportionlineage1=proportionlineage1_by_subject{i};
        isolate_compartments=compartment_lists{k};
        p1=nan(1,7);
        sem1=nan(1,7);
        p1(1)=mean(proportionlineage1(isolate_compartments==eta_site_number));
        sem1(1)=std(proportionlineage1(isolate_compartments==eta_site_number))/sqrt(sum(isolate_compartments==eta_site_number));
        fprintf(fid_F3,['P' num2str(k) ',' 'Tracheal aspirate,' num2str(p1(1)) ',' num2str(sem1(1)) '\n']);
        for c=1:6
            if sum(isolate_compartments==c)>0
                p1(c+1)=mean(proportionlineage1(isolate_compartments==c));
                sem1(c+1)=std(proportionlineage1(isolate_compartments==c))/sqrt(sum(isolate_compartments==c));
                fprintf(fid_F3,['P' num2str(k) ',' sites{c} ',' num2str(p1(c+1)) ',' num2str(sem1(c+1)) '\n']);
            end
        end
        if mean(p1)>.5
            p1=1-p1;
        end
        plot([0 8],[mean(p1(~isnan(p1))) mean(p1(~isnan(p1)))],'k-')
        fprintf(fid_F3,['P' num2str(k) ',' 'Mean across samples,' num2str(mean(p1(~isnan(p1)))) '\n']);
        bar([p1; 1-p1]','stacked','Edgecolor','None')
        colormap(cell2mat(cellfun(@rgb,{'DimGray','LightGray'},'UniformOutput', false)'))
        if i > 1
            plot([1:numel(p1); 1:numel(p1)], [p1-sem1; p1+sem1], 'k-','LineWidth',2)
        end
        
        set(gca,'xtick',[])
        set(gca,'ytick',0:.25:1)
        set(gca,'yticklabel',{'0', '','.5', '','1'});
        i=i+1;
        xlim([0 8])
    end
end



%% average coverage stats

number_subjects_for_which_samples_are_diverse_scrapes=12;

avgcovscrapes=sum(average_coverage_subject(1:number_subjects_for_which_samples_are_diverse_scrapes).*samples_per_subject(1:number_subjects_for_which_samples_are_diverse_scrapes))/sum(samples_per_subject(1:number_subjects_for_which_samples_are_diverse_scrapes));
avgcovcolonies=sum(average_coverage_subject(number_subjects_for_which_samples_are_diverse_scrapes+1:end).*samples_per_subject(number_subjects_for_which_samples_are_diverse_scrapes+1:end))/sum(samples_per_subject(number_subjects_for_which_samples_are_diverse_scrapes+1:end));

fprintf(1,['Avg cov scrapes ' num2str(avgcovscrapes) '\n']);
fprintf(1,['Avg cov colonies ' num2str(avgcovcolonies) '\n']);

%% compare dMRCA across groups
load([masterdir '/tools/prev_HIV_dx'])
% -1 indicates dx of HIV upon death; 0 indicates HIV negative; 1 indicates HIV positive

dMRCA_prev_HIV_dx=distance_to_subject_MRCA(prev_HIV_dx > 0 & ~ismember(subjects,subjects_with_multiple_strains)');
dMRCA_no_prev_HIV_dx=distance_to_subject_MRCA(prev_HIV_dx <1 & ~ismember(subjects,subjects_with_multiple_strains)');

p_value_different = ranksum(dMRCA_no_prev_HIV_dx,dMRCA_prev_HIV_dx);


%% Figure 2a & 2b


%dMRCA
figure; clf;
subplot(2,1,1)
hold on;
bins=[-1.65:.25:3.5];
plot_distance_to_subject_MRCA=distance_to_subject_MRCA;
plot_distance_to_subject_MRCA(plot_distance_to_subject_MRCA<.08)=.08;

hist(log10(plot_distance_to_subject_MRCA),bins);
xlabel('Average SNP distance to subject MRCA')
ylabel('Number of subjects')
h = findobj(gca,'Type','patch');
set(h,'FaceColor',rgb('Black'))
set(h,'EdgeColor',rgb('White'))
set(gca,'Ytick',0:5:20)
set(gca,'Ytick',0:5:20)
ticks=-1:.5:3.5;
set(gca,'Xtick',ticks(1:2:end))
set(gca,'Xticklabel',10.^ticks(1:2:end))
ylim([-.01 15])
xlim([-1.35 inf])

fid_F2=fopen([masterdir '/csvfiles_generated/figure2sourcedata.csv'],'w');
fprintf(fid_F2,'\n#Figure 2A\n#Also see Supplementary Table 2; Supplementary Table 4; and Supplementary Figure 3\n');
fprintf(fid_F2,'#Subject,Distance to subject MRCA\n');
for i=1:numel(plot_distance_to_subject_MRCA)
    fprintf(fid_F2,['P' num2str(i) ',' num2str(plot_distance_to_subject_MRCA(i)) '\n']);
end

% number of de novo mutations per subject
subplot(2,1,2); hold on;
[number_mutations_per_subject_n,subject_i]=find_count_duplicates([de_novo_muts_all.subject]);

number_mutations_per_subject=zeros(numel(subjects),1);
number_mutations_per_subject(subject_i)=number_mutations_per_subject_n;

hist(number_mutations_per_subject,0:1:63);
h = findobj(gca,'Type','patch');
set(h,'FaceColor',rgb('Black'))
set(h,'EdgeColor','white')
xlabel('Number of de novo mutations')
ylabel('Number of subjects')
set(gca,'Ytick',0:2:6)
xlim([-.75 66])
ylim([-.01 7])

fprintf(fid_F2,'\n\n#Figure 2B\n#Also see Supplementary Table 2; Supplementary Table 4; and Supplementary Figure 3\n');
fprintf(fid_F2,'#Subject,Number of de novo mutations per subject\n');
for i=1:numel(number_mutations_per_subject)
    fprintf(fid_F2,['P' num2str(i) ',' num2str(number_mutations_per_subject(i)) '\n']);
end


%% suppelemental figure about sampling

figure; clf;
subplot(2,2,1);hold on;
plot(samples_per_subject, number_mutations_per_subject,'ko','MarkerFaceColor',rgb('DarkGray'),'MarkerSize', 6)
plot(samples_per_subject(prev_HIV_dx==-1), number_mutations_per_subject(prev_HIV_dx==-1),'ko','MarkerFaceColor',rgb('Blue'),'MarkerSize', 6)
plot(samples_per_subject(prev_HIV_dx==0), number_mutations_per_subject(prev_HIV_dx==0),'ko','MarkerFaceColor',rgb('Red'),'MarkerSize', 6)
coeff=polyfit(samples_per_subject',number_mutations_per_subject',1);
plot(0:150,(0:150)*coeff(1)+coeff(2),'k-')
xlim([0 150])
set(gca,'Ytick',0:15:75)
a=corrcoef(samples_per_subject',number_mutations_per_subject');
disp(a(1,2)^2);

subplot(2,2,2);hold on;
plot(number_specimens, number_mutations_per_subject,'ko','MarkerFaceColor',rgb('DarkGray'),'MarkerSize', 6)
plot(number_specimens(prev_HIV_dx==-1), number_mutations_per_subject(prev_HIV_dx==-1),'ko','MarkerFaceColor',rgb('Blue'),'MarkerSize', 6)
plot(number_specimens(prev_HIV_dx==0), number_mutations_per_subject(prev_HIV_dx==0),'ko','MarkerFaceColor',rgb('Red'),'MarkerSize', 6)
coeff=polyfit(number_specimens',number_mutations_per_subject',1);
plot(0:11,(0:11)*coeff(1)+coeff(2),'k-')
xlim([2 11])
set(gca,'Ytick',0:15:75)
a=corrcoef(number_specimens',number_mutations_per_subject');
disp(a(1,2)^2);

simpleinf= ~ismember(subjects,subjects_with_multiple_strains)';

subplot(2,2,3);hold on;
plot(samples_per_subject(simpleinf), distance_to_subject_MRCA(simpleinf),'ko','MarkerFaceColor',rgb('Gray'),'MarkerSize', 6)
plot(samples_per_subject(prev_HIV_dx==-1 & simpleinf), distance_to_subject_MRCA(prev_HIV_dx==-1 & simpleinf),'ko','MarkerFaceColor',rgb('Blue'),'MarkerSize', 6)
plot(samples_per_subject(prev_HIV_dx==0 & simpleinf), distance_to_subject_MRCA(prev_HIV_dx==0 & simpleinf),'ko','MarkerFaceColor',rgb('Red'),'MarkerSize', 6)
coeff=polyfit(samples_per_subject(simpleinf)',distance_to_subject_MRCA(simpleinf)',1);
plot(0:150,(0:150)*coeff(1)+coeff(2),'k-')
xlim([0 150])
a=corrcoef(samples_per_subject(simpleinf)',distance_to_subject_MRCA(simpleinf)');
disp(a(1,2)^2);

subplot(2,2,4);hold on;
plot(number_specimens(simpleinf), distance_to_subject_MRCA(simpleinf),'ko','MarkerFaceColor',rgb('Gray'),'MarkerSize', 6)
plot(number_specimens(prev_HIV_dx==-1 & simpleinf), distance_to_subject_MRCA(prev_HIV_dx==-1& simpleinf),'ko','MarkerFaceColor',rgb('Blue'),'MarkerSize', 6)
plot(number_specimens(prev_HIV_dx==0 & simpleinf), distance_to_subject_MRCA(prev_HIV_dx==0 & simpleinf),'ko','MarkerFaceColor',rgb('Red'),'MarkerSize', 6)
coeff=polyfit(number_specimens(simpleinf)',distance_to_subject_MRCA(simpleinf)',1);
plot(0:11,(0:11)*coeff(1)+coeff(2),'k-')
xlim([2 11])
a=corrcoef(number_specimens(simpleinf)',distance_to_subject_MRCA(simpleinf)');
disp(a(1,2)^2);


%% proportion minor lineage

for i=1:numel(proportionlineage1_by_subject)
    x=1-mean(proportionlineage1_by_subject{i});
    if x > .6
        x=1-x;
    end
    fprintf(['Proportion Lineage 1: ' num2str(x) '\n']);
end


%% supplemental fig 2a legened


k=20;
figure; clf; hold on;

s=0:.1:1;
for i=1:numel(s)
    x=1;
    y=i;
    patch([x x x+1 x+1],[y-.5 y+.5 y+.5 y-.5],'r', 'FaceColor', (1- s(i))*[1 1 1],'EdgeColor', rgb('Black'))
end
set(gca,'Xtick',[])
set(gca,'Ytick',[])

%% colored block labels for supplemental figure 2c

%These orders copied in from .tree file
%P20
treeorder={'Liver-C5-7','Liver-C5-2','Liver-C4-1','Liver-C2-1','Liver-C1-4','Liver-A5-1','Liver-A4-1','Liver-A3-1','Liver-A1-1','Lung6-C4-1','Lung6-B5-1','Lung6-B4-1','Lung6-B2-4','Lung6-B1-1','Lung6-A3-1','Lung5-C4-1','Lung5-C2-3','Lung5-C1-1','Lung5-B5-1','Lung5-A5-1','Lung5-A1-5','Lung4-A4-7','Lung4-A4-2','Lung3-C3-1','Lung3-C2-3','Lung3-B5-4','Liver-C1-5','Lung6-B3-1','Lung6-B2-5','Lung5-C5-1','Lung5-A2-9','Lung5-A1-4','Lung4-C4-1','Lung4-C3-1','Lung4-A5-5','Lung4-A3-3','Lung4-C2-1','Lung4-A5-4','Lung4-A3-6','Lung4-A2-1','Lung4-A1-1','Lung3-C5-2','Lung3-C2-6','Lung3-C1-4','Lung3-B5-3','Lung3-B1-4','Liver-B3-6','Liver-B1-1','Lung5-C3-1','Lung5-C2-6','Lung3-C5-7','Lung3-C4-1','Lung3-C1-5','Lung3-B5-2','Lung3-B4-1','Lung3-B2-1','Lung3-B1-5','Lung3-A1-1','Lung2-C5-3','Lung2-C4-3','Lung2-C1-1','Lung2-B6-7','Lung2-B5-6','Lung2-B2-4','Lung2-A4-2','Lung2-A3-3','Lung2-A2-7','Lung2-A1-7','Lung2-C5-6','Lung2-C4-6','Lung2-B6-2','Lung2-B5-3','Lung2-B4-9','Lung2-B3-9','Lung2-B2-5','Lung2-B1-8','Lung2-A4-7','Lung2-A3-6','Lung2-A2-2','Lung2-A1-2','Lung1-B2-1','Lung1-A4-1','Lung1-A3-1','Lung1-A2-1','Lung1-A1-1','Eta-C4-100','Eta-B5-100','Eta-B4-100','Eta-B1-100','Eta-A5-100','Eta-A4-100','Eta-A3-100','Lung2-A5-9','Eta-B2-100','Eta-A1-90a'};
%P15
treeorder={'Lung6-C2-9','Lung5-C5-7','Lung5-C4-7','Lung5-C3-9','Lung5-C2-6','Lung5-B5-9','Lung5-B4-9','Lung5-B3-7','Lung5-B3-2','Lung5-A4-9','Lung5-A3-3','Lung5-A2-9','Lung4-C5-9','Lung4-C4-9','Lung4-C3-9','Lung4-C2-9','Lung4-C1-9','Lung4-B5-5','Lung4-B5-3','Lung4-B3-9','Lung4-B2-9','Lung4-B1-8','Lung4-A5-9','Lung4-A4-9','Lung4-A3-9','Lung4-A2-3','Lung4-B4-9','Lung4-A2-6','Lung3-C2-8','Lung3-B1-9','Lung3-A5-8','Lung3-A2-9','Lung2-C4-2','Lung2-B5-9','Lung2-B4-7','Lung2-B2-3','Lung5-A3-4','Lung2-A5-2','Lung2-A4-3','Lung2-A3-8','Lung2-A1-5','Eta-C4-98a','Eta-C3-99a','Eta-C1-99a','Eta-B5-98a','Eta-B4-98a','Eta-B2-98a','Eta-B1-98a','Eta-A4-97a','Eta-A3-98a','Eta-A1-97a','Liver-C4-9','Liver-C3-8','Liver-C2-5','Liver-C2-4','Liver-C1-8','Liver-A2-9','Liver-A1-9','Liver-C5-9','Lung6-C5-9','Lung6-C3-9','Lung6-B5-9','Lung6-B4-9','Lung6-B3-9','Lung6-B2-9','Lung6-B1-9','Lung6-A5-8','Lung6-A4-9','Lung6-A3-9','Lung6-C4-9','Lung6-C1-9','Lung6-A1-8','Lung5-C5-2','Lung5-C4-2','Lung5-C2-3','Lung5-A1-9','Lung4-A1-5','Lung4-A1-3','Lung3-C1-9','Lung3-A6-7','Lung3-A4-9','Lung3-A3-6','Lung3-A6-2','Lung3-A3-3','Lung2-C4-7','Lung2-C3-8','Lung2-C2-8','Lung2-C1-6','Lung2-C1-2','Lung2-B4-2','Lung2-B3-9','Lung2-B2-6','Lung2-B1-8','Lung2-A5-6','Lung2-A4-6','Lung2-A2-9','Lung2-A1-4','Lung1-C5-9','Lung1-C4-9','Lung1-C3-9','Lung1-C2-9','Lung1-C1-9','Lung1-B5-9','Lung1-B4-9','Lung1-B3-9','Lung1-B2-9','Lung1-B1-9','Lung1-A5-9','Lung1-A4-9','Lung1-A3-9','Lung1-A2-9','Lung1-A1-9','Eta-C5-98a','Eta-C2-98a','Eta-B3-98a','Eta-A2-98a'};
sitesk=zeros(1,numel(treeorder));
figure; clf; hold on;
for s=1:numel(treeorder)
    x=strsplit(treeorder{s},'-');
    plot([0 1],[s s], '-', 'Color', rgb(site_colors(strcmp(x{1},sites))),'LineWidth',3);
end
ylim([-2 numel(treeorder)+1])


figure; clf; hold on;
for s=1:numel(sites)
    plot([0 1],[s s], '-', 'Color', rgb(site_colors(s)),'LineWidth',10);
end


%% prepare for dNdS analysis

dNdSfilename=[REFGENOMEFOLDER '/dNdStools_' num2str(promotersize) '.mat'];
if exist(dNdSfilename,'file')
    load(dNdSfilename)
else
    div_mutation_type_probability_matrix(REFGENOMEFOLDER, promotersize);
    save(dNdSfilename, 'genomewide_possibilities', 'cds_possibilities')
end

%%  load in intersubject muts for compairsion for dNdS

cd([masterdir '/subject_folders/intersubject'])
load('de_novo_muts')

interpatient_genomic_positions=[de_novo_muts(:).pos];

interpatient_snps=interpatient_genomic_positions(~ismember(interpatient_genomic_positions,genomic_positions_all));
interpatient_types=[de_novo_muts(:).type];


Ni=sum(interpatient_types=='N');
Si=sum(interpatient_types=='S');

intersubject_typecounts=zeros(4,1);
intersubject_mutationmatrix=zeros(4,4);

for i=1:numel(interpatient_snps);
    
    anc=find(NTs==de_novo_muts(i).anc);
    new=de_novo_muts(i).nts;
    
    if find(new==NTs(anc) & numel(new)>1)
        new(new==NTs(anc))=[];
        new(new=='N')=[];
        new=find('ATCG'==new(1));
        
        intersubject_mutationmatrix(anc,new)=intersubject_mutationmatrix(anc,new)+1;
        
        if de_novo_muts(i).type=='N'
            intersubject_typecounts(1)=intersubject_typecounts(1)+1;
        elseif de_novo_muts(i).type=='S'
            intersubject_typecounts(2)=intersubject_typecounts(2)+1;
        elseif de_novo_muts(i).type=='P'
            intersubject_typecounts(3)=intersubject_typecounts(3)+1;
        elseif de_novo_muts(i).type=='I'
            intersubject_typecounts(4)=intersubject_typecounts(4)+1;
        else
            error('unrecognized type')
        end
    end
end





%% dndS  different number of mutations


categories={find(~cellfun(@isempty,{de_novo_muts_all.locustag}))} ;

mutationmatrix=zeros(4,4,numel(categories)); %last dimension is context
typecounts=zeros(numel(categories),4); %NSPI %first dimension is context


for context=1:numel(categories)
    
    muts_in_category=categories{context};
    
    for j=1:numel(muts_in_category);
        
        i=muts_in_category(j);
        
        anc=find(NTs==de_novo_muts_all(i).anc);
        new=de_novo_muts_all(i).nts;
        
        if find(new==NTs(anc))
            new(new==NTs(anc))=[];
            new(new=='N')=[];
            new=find('ATCG'==new);
            
            
            mutationmatrix(anc,new,context)=mutationmatrix(anc,new,context)+1;
            
            if de_novo_muts_all(i).type=='N'
                typecounts(context,1)=typecounts(context,1)+1;
            elseif de_novo_muts_all(i).type=='S'
                typecounts(context,2)=typecounts(context,2)+1;
            elseif de_novo_muts_all(i).type=='P'
                typecounts(context,3)=typecounts(context,3)+1;
            elseif de_novo_muts_all(i).type=='I'
                typecounts(context,4)=typecounts(context,4)+1;
            else
                error('unrecognized type')
            end
        end
    end
end

typecounts=[intersubject_typecounts'; typecounts];
mutationmatrix=cat(3,intersubject_mutationmatrix, mutationmatrix);

mutO = div_matrix2_6types(mutationmatrix);

percentN_matrix=(genomewide_possibilities(:,:,1)./(genomewide_possibilities(:,:,1)+genomewide_possibilities(:,:,2)));
percentN_types=div_matrix2_6types(percentN_matrix)./2;

No=typecounts(:,1);
So=typecounts(:,2);

expectedN=mutO'*percentN_types;
expectedS=sum(mutO)'-expectedN;

NSe=expectedN./expectedS;
NSo=typecounts(:,1)./typecounts(:,2);

mycolors=([rgb('green'); rgb('purple'); rgb('yellow'); rgb('red'); rgb('blue'); rgb('black'); rgb('gray')]);

%Indidiviudal
[upperCI,lowerCI]=binomialCIdNdS(No,So,NSe(1));
figure; clf; hold on;
for i=1:numel(NSo)
    h=bar(i,NSo(i)/NSe(i));
    set(h,'FaceColor',mycolors(i,:));
    %    h(1).BaseValue = 1;
end
plot([1:numel(NSo);1:numel(NSo)],[lowerCI';upperCI'],'k');
ylim([.4 2])
set(gca,'Ytick',[0.25 .5 1 2])
set(gca,'Yticklabel',{'0.25', '.5' '1' '2' '4'})
set(gca,'YScale','log')
set(gca,'Xtick',1:numel(NSo))
set(gca,'Xticklabel',{'Intersubject SNPs', 'Intrasubject SNPS', 'Genes mutated more than once', 'Genes mutated three times', 'Genes mutated twice within a subject'})

p_purifying_selection=binocdf(typecounts(:,1),typecounts(:,1)+typecounts(:,2),expectedN./(expectedN+expectedS));



figure; clf; hold on;
h=bar(1,NSo(2)/NSe(2));
set(h,'FaceColor',rgb('gray'));
plot([1 1],[lowerCI(2) upperCI(2)],'k');
h=bar(2,NSo(1)/NSe(1));
set(h,'FaceColor',rgb('white'));
plot([2 2],[lowerCI(1) upperCI(1)],'k');
ylim([.6 1.5])
xlim([0 3])
set(gca,'Ytick',[0.25 .5 .75 1 1.5 2 4])
set(gca,'Xtick',[4])
set(gca,'YScale','log')




%% show mutation spectra


mutationmatrix_all=zeros(4,4);

for i=1:numel(de_novo_muts_placed)
    anc=find(NTs==de_novo_muts_placed(i).anc);
    new=de_novo_muts_placed(i).nts;
    if find(new==NTs(anc))
        new(new==NTs(anc))=[];
        new(new=='N')=[];
        new=find('ATCG'==new);
        mutationmatrix_all(anc,new)=mutationmatrix_all(anc,new)+1;
    end
end
mutO_interpateint=div_matrix2_6types(intersubject_mutationmatrix)/sum(intersubject_mutationmatrix(:));
mutO_intrasubject=div_matrix2_6types(mutationmatrix_all)/sum(mutationmatrix_all(:));




%build expecation based on %normalize for GC content
ATCGpercentongenome=sum(sum(genomewide_possibilities,3),2)/sum(genomewide_possibilities(:));
GCcontent=ATCGpercentongenome(3)+ATCGpercentongenome(4);
ATcontent=1-GCcontent;

%expectation based on transitions being twice as likely as transversions
expectation=[ATcontent/4 ATcontent/4 GCcontent/4 GCcontent/4 ATcontent/2 GCcontent/2];


figure; clf; hold on;
bar([mutO_intrasubject([1 2 4 5 3 6]) mutO_interpateint([1 2 4 5 3 6])], 'grouped'); colormap([rgb('Gray'); rgb('White');]);
plot([(1:6)-.25; (1:6)+.3],[expectation; expectation], 'k--', 'LineWidth', 2);
ylabel('Percent of SNPs')
set(gca,'Ytick',0:.1:.5)
set(gca,'Xtick',1:6)
set(gca,'Xticklabel',{'AT->TA', 'AT->CG', 'GC->CG', 'GC->TA', 'AT->GC', 'GC->AT'})
%legend({'Intrasubject SNPs', 'Intersubject SNPs'})
ylim([0 .5])


a=div_matrix2_6types(intersubject_mutationmatrix);
b=div_matrix2_6types(mutationmatrix_all);

pvalue_for_intra_vs_inter_spectrum=fexact([b(5) b(6); a(5) a(6);]);


% show probability amino acid changing
figure; clf; hold on;
bar(percentN_types); colormap('Gray');
set(gca,'Xtick',1:6)
set(gca,'Xticklabel',{}) %{'AT->TA', 'AT->CG', 'GC->CG', 'GC->TA', 'AT->GC', 'GC->AT'})
set(gca,'Ytick',.5:.1:1)
ylim([.5 .9])
ylabel('Probability amino acid changing')




%% show how diveristy within lung is distributed, include eta -- Figure 3C

fprintf(fid_F3,'\n\n#Figure 3C\n#Also see Supplementary Table 4 for more information on each mutated position and distribution across samples\n');

subjects_to_plots=[20 41];
colorlist_t={'Black', 'LightGrey', 'LightCoral', 'CornflowerBlue', 'Cyan', 'Peru', 'Brown', 'Violet', 'Orange', 'MediumBlue', 'MediumPurple', 'LightPink', 'ForestGreen', 'PaleGreen','MediumSeaGreen','Black', 'LightGrey', 'CornflowerBlue', 'LightCoral', 'Cyan', 'Peru', 'Brown', 'Violet', 'Orange', 'MediumBlue', 'MediumPurple', 'LightPink', 'ForestGreen', 'PaleGreen','MediumSeaGreen'};


for k=subjects_to_plots
    
    fprintf(fid_F3,['\n#Subject P' num2str(k) '\n']);
    fprintf(fid_F3,'#Genomic positions mutated in genotype(semicolon seperated), Frequency in tracheal aspirate,Frequency in Lung1, Frequency in Lung2, Frequency in Lung3, Frequency in Lung4, Frequency in Lung5, Frequency in Lung6\n');
    
    isolate_compartments=compartment_lists{k};
    genotypes=genotypematrices{k};
    de_novo_muts_pt=de_novo_muts_placed([de_novo_muts_placed(:).subject]==k);
    
    if ~isempty(genotypes) & sum(sum(genotypes_freq_sites_matrices{k}(:,1:eta_site_number))>0)>1
        
        %normalize
        sitefreqs=genotypes_freq_sites_matrices{k};
        sitefreqs(:,(eta_site_number+1):end)=0;
        sitefreqs=sitefreqs./repmat(sum(sitefreqs),size(sitefreqs,1),1);
        
        %remove low frequency strains from this plotting
        genotypestokeep=sum(sitefreqs(1:end-1,:)>= min_freq_within_site,2) >0;
        genotypestokeep=[genotypestokeep; true];
        plotsitefreqs = sitefreqs(genotypestokeep,:) ./ repmat(sum(sitefreqs(genotypestokeep,:)),sum(genotypestokeep),1);
        plotgenotypes=[genotypes(genotypestokeep(1:end-1),:); zeros(1,size(genotypes,2))]; %keep ancestor
        
        
        if ~isempty(genotypes)
            
            %reorder genotypes by number of mutations
            forsorting=sum(plotgenotypes,2)+400*plotgenotypes(:,end);
            [~,sortedorder]=sort(forsorting,'ascend');
            neworder=sortedorder;
            
            if k==41
                neworder=sortedorder([1 3 4 6 8 10 11 7 9 5 12 13 14 2]);
                colorlist_t={'Black', 'Red', 'Salmon','Orange','DarkKhaki','DarkKhaki','DarkKhaki','Gold', ...
                    'Khaki','Lavender','PaleGreen',... %'GreenYellow'...
                    'Indigo' ... %'Turquoise', 'DeepSkyBlue',
                    'MidnightBlue','Sienna'};
            elseif k==20
                colorlist_t={'Black', 'FireBrick', 'Salmon', 'Maroon','Orange','Chocolate','PaleGreen','MediumSeaGreen'}; ...
                    %'Turquoise', 'MediumBlue'}; %, 'MediumBlue', 'MediumPurple', 'LightPink', 'PaleGreen','MediumSeaGreen','Black', 'LightGrey', 'CornflowerBlue', 'LightCoral', 'Cyan', 'Peru', 'Brown', 'Violet', 'Orange', 'MediumBlue', 'MediumPurple', 'LightPink', 'ForestGreen', 'PaleGreen','MediumSeaGreen'};
            else
                neworder=sortedorder;
            end
            
            %implement new order
            plotgenotypes=plotgenotypes(neworder,:);
            plotsitefreqs=plotsitefreqs(neworder,:);
            
            %tree
            nts_for_tree=char(plotgenotypes);
            nts_for_tree(plotgenotypes==0)='A';
            nts_for_tree(plotgenotypes==1)='T';
            
            TreeSampleNames2={};
            TreeColors={};
            for i=1:sum(genotypestokeep)
                TreeSampleNames2{end+1}=['GT' num2str(i)];
                TreeColors{end+1}='16777216'; % all same for now
            end
            
            cd([masterdir '/trees/temptreefolder'])
            treefilename=generate_parsimony_tree(nts_for_tree', TreeSampleNames2, ['P' num2str(k)]);
            color_tree(treefilename, [masterdir '/trees/genotypetrees_lung/P' num2str(k) '.tree'], TreeSampleNames2, TreeColors)
            cd(masterdir)
            
            sitesk=find(sum(plotsitefreqs)>0);
            
            %put endotracheal aspirate first
            if ismember(eta_site_number,sitesk)
                sitesk(sitesk==eta_site_number)=[];
                sitesk=[eta_site_number sitesk];
            end
            
            %plot pie charts as stack bars
            figure; clf; hold on;
            plotsitefreqs=plotsitefreqs(:,sitesk);
            bar(plotsitefreqs(:,sum(plotsitefreqs)>0)','stacked','Edgecolor','None')
            colormap(cell2mat(cellfun(@rgb,colorlist_t,'UniformOutput', false)'))
            set(gca,'xtick',[])
            set(gca,'ytick',0:.25:1)
            set(gca,'yticklabel',{'0','','.5','','1'})
            xlim([0 8])
            
            figure; clf; hold on;
            
            for i=1:numel(colorlist_t)
                plot(5,50-i,'s','MarkerFaceColor',rgb(colorlist_t{i}),'MarkerEdgeColor','None','MarkerSize',15)
                text(5.1,50-i,num2str(i))
            end
            ylim([50-i-5 50])
            
            for i=1:size(plotsitefreqs,1)
                muts_genotype=find(plotgenotypes(i,:));
                if isempty(muts_genotype)
                    fprintf(fid_F3,'Anc');
                else
                    for j=1:numel(muts_genotype);
                        fprintf(fid_F3,[num2str(de_novo_muts_pt(muts_genotype(j)).pos) ';']);
                    end
                end
                fprintf(fid_F3,',');
                for j=1:numel(sitesk)
                    fprintf(fid_F3,[num2str(plotsitefreqs(i,j)) ',']);
                end
                fprintf(fid_F3,'\n');
            end
            
            
        end
    end
end



%% pairwise distance analysis, don't include eta when thinking about intraorgan

intrasitemeans_means=nan(numel(subjects),1);
intersitemeans_means=nan(numel(subjects),1);
interorganmeans_means=nan(numel(subjects),1);

intrasitemeans_SEM=nan(numel(subjects),1);
intersitemeans_SEM=nan(numel(subjects),1);
interorganmeans_SEM=nan(numel(subjects),1);

inter_site_pair_means=nan(numel(subjects),numel(sites),numel(sites));
intra_site_means_per_site=nan(numel(subjects),numel(sites));
pvalues=nan(numel(subjects),4);



fid_F4=fopen([masterdir '/csvfiles_generated/figure4sourcedata.csv'],'w');
fprintf(fid_F3,'\n\n#Figure 3B\n#Also see Supplementary Table 4\n');
fprintf(fid_F3,'#Only subjects with multiple samples per site considered\n');
fprintf(fid_F3,'#Subject,Intrasite SNP distance (mean across sites),Intrasite SNP distance (mean each site - semicolon delimited),Intersite SNP distance (mean across pairs),Intersite SNP distance (mean each pair - semicolon delimited) \n');

fprintf(fid_F4,'\n#Figure 4B\n#Also see Supplementary Table 4 and Supplementary Figure 6\n');
fprintf(fid_F4,'#Only subjects with multiple samples per site considered\n');
fprintf(fid_F4,'#Subject,Intersite SNP distance (mean across pairs),Intersite SNP distance (mean each pair - semicolon delimited), Interorgan SNP distance (mean across pairs),Interorgan SNP distance (mean each pair - semicolon delimited) \n');


for k=1:numel(subjects)
    
    
    
    %initialize
    site_list=compartment_lists{k};
    [number_samples_site,sites_found]=find_count_duplicates(site_list);%only include sites sampled more than once
    sites_pt=sites_found(number_samples_site>1);
    lungsites=sites_pt(sites_pt<eta_site_number);
    ptmatrix=observed_matrices_only_de_novo_mutations{k}(:,ismember(site_list,sites_pt));
    site_list=site_list(ismember(site_list,sites_pt));
    
    %calculate distance for intrasite pairs
    intrasitepairs=[];
    instrasitemean=0;
    for c=1:numel(sites_pt)
        samples=find(site_list==sites_pt(c));
        instrasitemean=0;
        for i=1:numel(samples)
            p1=ptmatrix(:,samples(i));
            for j=(i+1):numel(samples)
                p2=ptmatrix(:,samples(j));
                pairwisedistance=sum(p1)+sum(p2)-sum(2*(p1.*p2));
                instrasitemean=instrasitemean+pairwisedistance;
                if sites_pt(c) < eta_site_number
                    intrasitepairs(end+1)=pairwisedistance;
                end
            end
        end
        intra_site_means_per_site(k,sites_pt(c))=instrasitemean/(.5*numel(samples)*(numel(samples)-1));
    end
    intrasitemeans=intra_site_means_per_site(k,~isnan(intra_site_means_per_site(k,1:(eta_site_number-1))));
    
    
    %intersite pairs
    intersitepairs=[];
    interorganpairs=[];
    for c1=1:numel(sites_pt)
        samples1=find(site_list==sites_pt(c1));
        for c2=(c1+1):numel(sites_pt)
            samples2=find(site_list==sites_pt(c2));
            intersitemean=0;
            for i=1:numel(samples1)
                p1=ptmatrix(:,samples1(i));
                for j=1:numel(samples2)
                    p2=ptmatrix(:,samples2(j));
                    pairwisedistance=sum(p1)+sum(p2)-sum(2*(p1.*p2));
                    intersitemean=intersitemean+pairwisedistance;
                    if sites_pt(c1) < eta_site_number & sites_pt(c2) < eta_site_number  %only save lung sites here
                        intersitepairs(end+1)=pairwisedistance;
                    elseif sites_pt(c1) < eta_site_number & ...
                            (sites_pt(c2)==liver_site_number | sites_pt(c2)==spleen_site_number)
                        interorganpairs(end+1)=pairwisedistance;
                    end
                end
            end
            inter_site_pair_means(k,sites_pt(c1),sites_pt(c2))=intersitemean/(numel(samples1)*numel(samples2));
        end
    end
    
    
    intersitemeans=inter_site_pair_means(k,1:(eta_site_number-1),1:(eta_site_number-1));
    intersitemeans=intersitemeans(~isnan(intersitemeans));
    
    interorganmeans=inter_site_pair_means(k,1:(eta_site_number-1),[liver_site_number spleen_site_number]);
    interorganmeans=interorganmeans(~isnan(interorganmeans));
    
    
    
    %summary and stats
    % 1 & 2 --Intrasite vs. Between site
    % 3 & 4 --Between site vs. Between organs
    if numel(intrasitepairs)>0 & numel(intersitepairs)>0
        [~,p1]=ttest2(intrasitepairs,intersitepairs, 'Tail', 'left');
        [~,p2]=ttest2(intrasitemeans,intersitemeans, 'Tail', 'left');
        pvalues(k,1:2)=[p1 p2];
    end
    if numel(interorganpairs)>0 & numel(intersitepairs)>0
        [~,p3]=ttest2(intersitepairs,interorganpairs, 'Tail', 'left');
        [~,p4]=ttest2(intersitemeans,interorganmeans, 'Tail', 'left');
        pvalues(k,3:4)=[p3 p4];
    end
    
    intrasitemeans_means(k)=mean(intrasitemeans);
    intersitemeans_means(k)=mean(intersitemeans);
    interorganmeans_means(k)=mean(interorganmeans);
    
    intrasitemeans_SEM(k)=std(intrasitemeans)/sqrt(numel(intrasitemeans));
    intersitemeans_SEM(k)=std(intersitemeans)/sqrt(numel(intersitemeans));
    interorganmeans_SEM(k)=std(interorganmeans)/sqrt(numel(intersitemeans));
    
    if ~isnan(mean(intersitemeans)) & ~isnan(mean(intrasitemeans))
        fprintf(fid_F3,['P' num2str(k) ',' num2str(mean(intrasitemeans)) ',']);
        for i=1:numel(intrasitemeans)
            fprintf(fid_F3,[num2str(intrasitemeans(i)) ';']);
        end
        fprintf(fid_F3,[',' num2str(mean(intersitemeans)) ',']);
        for i=1:numel(intersitemeans)
            fprintf(fid_F3,[num2str(intersitemeans(i)) ';']);
        end
        fprintf(fid_F3,'\n');
    end
    
    if ~isnan(mean(interorganmeans))
        fprintf(fid_F4,['P' num2str(k) ',' num2str(mean(intersitemeans)) ',']);
        for i=1:numel(intersitemeans)
            fprintf(fid_F4,[num2str(intersitemeans(i)) ';']);
        end
        fprintf(fid_F4,[',' num2str(mean(interorganmeans)) ',']);
        for i=1:numel(interorganmeans)
            fprintf(fid_F4,[num2str(interorganmeans(i)) ';']);
        end
        fprintf(fid_F4,'\n');
    end
end



%% diversity by location of lung specimen


fprintf(fid_F4,'\n\n#Figure 4A\n#Also see Supplementary Figure 6 and Supplementary Table 4\n');
fprintf(fid_F4,'#Subject,Right lung intralung SNP distance (mean across site pairs),Right lung intralung SNP distance (mean each site pair - semicolon delimited), Left lung intralung SNP distance (mean across site pairs),Left lung intralung SNP distance (mean each site pair - semicolon delimited),Interlung SNP distance (mean across site pairs),Interlung SNP distance (mean each site pair - semicolon delimited)\n');

%site_pairs=find(sum(inter_site_pair_means)>0);
inter_site_pair_means_normalized=inter_site_pair_means;
intra_site_means_per_site_normalized=intra_site_means_per_site;

both_left_lung=zeros(numel(sites),numel(sites)); both_left_lung(1:3,1:3)=1;
both_right_lung= zeros(numel(sites),numel(sites)); both_right_lung(4:6,4:6)=1;
inter_lung_pairs= zeros(numel(sites),numel(sites)); inter_lung_pairs(1:6,1:6)=1; inter_lung_pairs(both_left_lung | both_right_lung)=0;
eta_site_pair= zeros(numel(sites),numel(sites)); eta_site_pair(1:6,8)=1;
any_lung=both_left_lung | both_right_lung | inter_lung_pairs | eta_site_pair ;

left_lung_mean=zeros(size(subjects));
right_lung_mean=zeros(size(subjects));
inter_lung_mean=zeros(size(subjects));


%normalize
for k=1:numel(subjects)
    
    intra_site_means_per_site_this_pt=intra_site_means_per_site(k,:);
    % only normalize if there are at least 3 lung sites from this subject
    if sum(~isnan(intra_site_means_per_site_this_pt))>2
        subject_average=mean(intra_site_means_per_site_this_pt(~isnan(intra_site_means_per_site_this_pt)));
        intra_site_means_per_site_normalized(k,:)=intra_site_means_per_site(k,:)/subject_average;
    else
        intra_site_means_per_site_normalized(k,:)=nan;
    end
    
    
    inter_site_pair_means_this_pt=squeeze(inter_site_pair_means(k,:,:));
    if sum(~isnan(intra_site_means_per_site_this_pt))>2
        % only normalize if there are at least 3 lung sites from this subject
        subject_average=mean(inter_site_pair_means_this_pt(any_lung & ~isnan(inter_site_pair_means_this_pt)));
        inter_site_pair_means_normalized(k,:,:)=inter_site_pair_means(k,:,:)/subject_average;
        
        left_lung_pairs_distance=inter_site_pair_means(k,both_left_lung & ~isnan(inter_site_pair_means_this_pt));
        right_lung_pairs_distance=inter_site_pair_means(k,both_right_lung & ~isnan(inter_site_pair_means_this_pt));
        inter_lung_pairs_distance=inter_site_pair_means(k,inter_lung_pairs & ~isnan(inter_site_pair_means_this_pt));
        
        left_lung_mean(k)=mean(left_lung_pairs_distance);
        right_lung_mean(k)=mean(right_lung_pairs_distance);
        inter_lung_mean(k)=mean(inter_lung_pairs_distance);
        
        
        if numel(right_lung_pairs_distance) > 0
            fprintf(fid_F4,['P' num2str(k) ',' num2str(mean(right_lung_pairs_distance)) ',']);
        else
            fprintf(fid_F4,['P' num2str(k) ',No pairs,']);
        end
        for i=1:numel(right_lung_pairs_distance)
            if ~isnan(right_lung_pairs_distance(i))
                fprintf(fid_F4,[num2str(right_lung_pairs_distance(i)) ';']);
            end
        end
        if numel(left_lung_pairs_distance) > 0
            fprintf(fid_F4,[',' num2str(mean(left_lung_pairs_distance)) ',']);
        else
            fprintf(fid_F4,',No pairs,');
        end
        for i=1:numel(left_lung_pairs_distance)
            if ~isnan(left_lung_pairs_distance(i))
                fprintf(fid_F4,[num2str(left_lung_pairs_distance(i)) ';']);
            end
        end
        if numel(inter_lung_pairs_distance) > 0
            fprintf(fid_F4,[',' num2str(mean(inter_lung_pairs_distance)) ',']);
        else
            fprintf(fid_F4,',No pairs,');
        end
        for i=1:numel(inter_lung_pairs_distance)
            if ~isnan(inter_lung_pairs_distance(i))
                fprintf(fid_F4,[num2str(inter_lung_pairs_distance(i)) ';']);
            end
        end
        fprintf(fid_F4,'\n');
    else
        inter_site_pair_means_normalized(k,:)=nan;
        left_lung_mean(k)=nan;
        right_lung_mean(k)=nan;
        inter_lung_mean(k)=nan;
    end
    
    
    
    
end


%% Same lung vs different lungs


h=figure(30); clf; %figure 4a
subplot(1,2,1); hold on;

set(h,'Position',[0 0 300 300]);
xlabel('Within same lung -- mean pairwise distance by subject')
ylabel('Different lungs -- mean pairwise distance by subject')
hold on; plot(0:.05:7,0:.05:7, 'k');
xlim([0 7])
ylim([0 7])
set(gca,'Xtick',0:1:7)
set(gca,'Ytick',0:1:7)
for i=1:numel(intersitemeans_means)
    if ~isnan(inter_lung_mean(i)) & ~isnan(right_lung_mean(i)) & ~isnan(left_lung_mean(i))
        plot([.01+left_lung_mean(i) .01+right_lung_mean(i)], [.01+inter_lung_mean(i) .01+inter_lung_mean(i)], '-','Color', rgb('Gray'),'LineWidth', .5)
        
        plot(.01+left_lung_mean(i), .01+inter_lung_mean(i),'ro', 'MarkerSize', 3, 'MarkerFaceColor', rgb('MediumBlue'),'MarkerEdgeColor', 'none');
        plot(.01+right_lung_mean(i), .01+inter_lung_mean(i),'bo', 'MarkerSize', 3, 'MarkerFaceColor', rgb('DodgerBlue'),'MarkerEdgeColor', 'none');
        
    end
end





%% localized diversity at site level

%Intrasite vs. Between site


h=figure(29);  %figure 3b
clf; subplot(1,2,1)
hold on;

set(h,'Position',[0 0 300 300]);
xlabel('within lung sites -- mean pairwise distance by subject')
ylabel('between lung sites -- mean pairwise distance by subject')
issignificant=~isnan(pvalues(:,1)) & pvalues(:,1)<.05;
plot(.001:.005:6,.001:.005:6, 'k');
xlim([0 5.5])
ylim([0 5.5])
set(gca,'Xtick',[0:1:6])
set(gca,'Ytick',[0:1:6])

% %error bars
for i=1:numel(intersitemeans_means)
    if ~isnan(intrasitemeans_means(i))
        plot([intrasitemeans_means(i)-intrasitemeans_SEM(i) intrasitemeans_means(i)+intrasitemeans_SEM(i)],[intersitemeans_means(i) intersitemeans_means(i)], '-', 'Color', rgb('Gray'),'LineWidth', .25);
        plot([intrasitemeans_means(i) intrasitemeans_means(i)],[intersitemeans_means(i)-intersitemeans_SEM(i) intersitemeans_means(i)+intersitemeans_SEM(i)], '-','Color', rgb('Gray'), 'LineWidth', .25);
    end
end
clickable_scatter_twocolor_oneoutput(intrasitemeans_means,intersitemeans_means, zeros(size(issignificant)))

%Mean x, SEM x, Mean x, SEM x,

%% color previous graph by hiv (comment out for graph in main text)

% plot(intrasitemeans_means(prev_HIV_dx>0),intersitemeans_means(prev_HIV_dx>0),'ko','MarkerFaceColor',rgb('Black'))
% plot(intrasitemeans_means(prev_HIV_dx<0),intersitemeans_means(prev_HIV_dx<0),'ko','MarkerFaceColor',rgb('Blue'))
% plot(intrasitemeans_means(prev_HIV_dx==0),intersitemeans_means(prev_HIV_dx==0),'ko','MarkerFaceColor',rgb('Red'))
%


%% Between site vs. Between organs

%figure 4b

figure(30); subplot(1,2,2); hold on;
set(h,'Position',[0 0 300 300]);
xlabel('between lung sites -- mean  pairwise distance by subject')
ylabel('between organs -- mean pairwise distance by subject')
hold on; plot(0:.5:6,0:.5:6, 'k');
issignificant=~isnan(pvalues(:,2)) & pvalues(:,2)<.001;
xlim([0 5.5])
ylim([0 5.5])
set(gca,'Xtick',0:1:5)
set(gca,'Ytick',0:1:5)
%figure; hold on;
for i=1:numel(intersitemeans_means)
    if ~isnan(interorganmeans_means(i))
        plot([intersitemeans_means(i)-intersitemeans_SEM(i) intersitemeans_means(i)+intersitemeans_SEM(i)],[interorganmeans_means(i) interorganmeans_means(i)], '-','Color', rgb('Gray'),'LineWidth', .25)
        plot([intersitemeans_means(i) intersitemeans_means(i)],[interorganmeans_means(i)-interorganmeans_SEM(i) interorganmeans_means(i)+interorganmeans_SEM(i)], '-', 'Color', rgb('Gray'),'LineWidth', .25)
    end
end
clickable_scatter_twocolor_oneoutput(intersitemeans_means,interorganmeans_means, zeros(size(issignificant)))

%% alpha diveristy across sites

lungsites_alpha_diversity=[];

figure; clf; hold on;
patch([3.5 6.5 6.5 3.5],[0 0 6 6],  rgb('LightGray'), 'EdgeColor', 'None');
%plot([0 12], [1 1] , 'k-', 'LineWidth', 1)

for p=1:spleen_site_number
    d=intra_site_means_per_site_normalized(:,p);
    d=d(~isnan(d));
    disp(numel(d))
    quantiles_d=quantile(d,[.25 .5 .75]) ;
    
    if p <eta_site_number
        lungsites_alpha_diversity=[lungsites_alpha_diversity; d];
        x=p;
    else
        x=p+1;
    end
    plot(x*ones(numel(d),1),d,'ko','MarkerFaceColor', rgb('Gray'), 'MarkerSize', 5)
    plot([x-.4 x+.4], [quantiles_d(1) quantiles_d(1)] , 'k-', 'LineWidth', 1)
    plot([x-.4 x+.4], [quantiles_d(2) quantiles_d(2)] , 'k-', 'LineWidth', 3)
    plot([x-.4 x+.4], [quantiles_d(3) quantiles_d(3)] , 'k-', 'LineWidth', 1)
    plot([x-.4 x-.4], [quantiles_d(1) quantiles_d(3)] , 'k-', 'LineWidth', 1)
    plot([x+.4 x+.4], [quantiles_d(1) quantiles_d(3)] , 'k-', 'LineWidth', 1)
    
end
set(gca,'Xtick',[1:eta_site_number])
set(gca,'Xticklabel',{});
set(gca,'Ytick',0:1:6)
xlabel ('Site')
ylabel('Alpha diveristy')


for p=1:eta_site_number
    d=intra_site_means_per_site_normalized(:,p);
    d=d(~isnan(d));
    disp(numel(d))
    plot(p*ones(numel(d),1),d,'ko','MarkerFaceColor', rgb('Gray'))
    quantiles_d=quantile(d,[.25 .5 .75]) ;
    plot([p-.4 p+.4], [quantiles_d(1) quantiles_d(1)] , 'k-', 'LineWidth', 1)
    plot([p-.4 p+.4], [quantiles_d(2) quantiles_d(2)] , 'k-', 'LineWidth', 3)
    plot([p-.4 p+.4], [quantiles_d(3) quantiles_d(3)] , 'k-', 'LineWidth', 1)
    plot([p-.4 p-.4], [quantiles_d(1) quantiles_d(3)] , 'k-', 'LineWidth', 1)
    plot([p+.4 p+.4], [quantiles_d(1) quantiles_d(3)] , 'k-', 'LineWidth', 1)
    if p <7
        lungsites_alpha_diversity=[lungsites_alpha_diversity; d];
    end
    
end


%% eta vs lung for diagnosis


fprintf(fid_F3,'\n\n#Figure 3D\n#Also see Supplementary Figure 4 and Supplementary Tables 2 and 4\n');
fprintf(fid_F3,'#Only subjects with at least one de novo mutation and a tracheal aspirate sample are included here\n');
fprintf(fid_F3,'#Subject,Number of de novo mutations observed,Number of de novo mutations observed in tracheal aspirate\n');

figure(29);
subplot(1,2,2);
hold on;
set(h,'Position',[0 0 300 300]);

%number of mutations in endotracheal aspirate
number_mutations_eta_per_subject=nan(size(subjects));
number_mutations_per_subject=nan(size(subjects));
for k=[20 41]; %1:numel(subjects)
    if sum(compartment_lists{k}==eta_site_number)>0 & sum(observed_matrices_only_de_novo_mutations{k}(:))>0
        etamatrix=observed_matrices_only_de_novo_mutations{k}(:,compartment_lists{k}==eta_site_number);
        number_mutations_eta_per_subject(k)= sum(sum(etamatrix,2)>.05);
        number_mutations_per_subject(k)=sum(sum(observed_matrices_only_de_novo_mutations{k},2)>.05);
        fprintf(fid_F3,['P' num2str(k) ',' num2str(number_mutations_per_subject(k)) ',' num2str(number_mutations_eta_per_subject(k)) '\n']);
        
    end
end

number_mutations_per_subject=number_mutations_per_subject(~isnan(number_mutations_eta_per_subject));
number_mutations_eta_per_subject=number_mutations_eta_per_subject(~isnan(number_mutations_eta_per_subject));

%jitter absolute
jittermax=.5;
number_mutations_eta_per_subject=max(number_mutations_eta_per_subject+(jittermax/2-jittermax*rand(1,numel(number_mutations_eta_per_subject))),0);
number_mutations_per_subject=max(number_mutations_per_subject+(jittermax/2-jittermax*rand(1,numel(number_mutations_per_subject))),0);


xlabel('de novo mutations in eta')
ylabel('de novo mutations in subject')
hold on; %plot(0:.5:4.5,0:.5:4.5, 'k');
plot(number_mutations_eta_per_subject,number_mutations_per_subject,'o','MarkerEdgeColor', 'none', 'MarkerFaceColor',rgb('Red'),'MarkerSize',4);
set(gca,'Xtick',0:10:60)
set(gca,'Ytick',0:10:60)
plot(.001:.005:45,.001:.005:45, 'k');
xlim([0 65])
ylim([0 65])


%% elaborate more on whats missing from endotracheal aspirate

%determine frequency of each mutation in each subject, averaging across sites
%Did the sputum proxy typically contain the most highly abundant variants?

freq_mutation_in_subject = [];
freq_mutation_in_eta = [];
found_in_eta = [];
number_specimens_in = [] ;
cutoff_within_sample=.05;
cutoff_within_specimen = .02;

for k=1:numel(subjects)
    if sum(compartment_lists{k}==eta_site_number)>0 & sum(observed_matrices_only_de_novo_mutations{k}(:))>0
        sites_pt=unique(compartment_lists{k});
        omatrix=observed_matrices_only_de_novo_mutations{k};
        omatrix_sites=zeros(length(sites_pt),size(omatrix,1));
        found_in_site=zeros(length(sites_pt),size(omatrix,1));
        
        for i=1:numel(sites_pt)
            c=sites_pt(i);
            omatrix_sites(i,:)=sum(omatrix(:,compartment_lists{k}==c),2)/sum(compartment_lists{k}==c);
            found_in_site(i,:)=sum(omatrix(:,compartment_lists{k}==c)>cutoff_within_sample,2)>0;
        end
        
        freq_mutation_in_subject = [freq_mutation_in_subject mean(omatrix_sites(sites_pt~=eta_site_number,:),1)];
        freq_mutation_in_eta = [freq_mutation_in_eta omatrix_sites(sites_pt==eta_site_number,:)];
        found_in_eta = [found_in_eta sum(omatrix(:,compartment_lists{k}==eta_site_number)>cutoff_within_sample,2)'];
        number_specimens_in = [number_specimens_in sum(found_in_site(sites_pt~=eta_site_number,:))];
    end
    
end

%remove positions only found in eta

found_in_eta=found_in_eta(number_specimens_in>0);
freq_mutation_in_eta=freq_mutation_in_eta(number_specimens_in>0);
freq_mutation_in_subject=freq_mutation_in_subject(number_specimens_in>0);
number_specimens_in=number_specimens_in(number_specimens_in>0);


%%  plot number of mutations observed and detected

cutoff=.05;

%probability detected by number of samples found in
figure; clf;
subplot(2,1,1);
hold on;

%calculate
bins=[0:.1:1];
num_eta=zeros(numel(bins)-1,1);
num_muts=zeros(numel(bins)-1,1);
for i=1:(numel(bins)-1)
    muts_in_bin= freq_mutation_in_subject>bins(i) & freq_mutation_in_subject <bins(i+1);
    num_eta(i) = sum(freq_mutation_in_eta(muts_in_bin)>cutoff);
    num_muts(i) = sum(muts_in_bin);
end
%plot
for i=2:numel(bins)-1
    patch([bins(i) bins(i) bins(i+1) bins(i+1)], [0 num_muts(i) num_muts(i) 0], rgb('White'));
    patch([bins(i) bins(i) bins(i+1) bins(i+1)], [0 num_eta(i) num_eta(i) 0], rgb('Gray'));
end
%make split y axis
maxy=max(num_eta(1:end));
maxyp=maxy+30;
ylim([0 maxyp+2])
patch([bins(1) bins(1) bins(2) bins(2)], [0 maxyp maxyp 0], rgb('White'));
patch([bins(1) bins(1) bins(2) bins(2)], [0 num_eta(1) num_eta(1) 0], rgb('Gray'));
set(gca,'Ytick',[0:5:20 maxyp])
set(gca,'Yticklabel',num2str([0:5:20 num_muts(1)]'))

%show cumulative probability with split bins for better visualization
prob_in_eta=zeros(numel(bins),1);
num_muts=zeros(numel(bins));
for i=1:numel(bins)
    muts_in_bin=  freq_mutation_in_subject > bins(i);
    prob_in_eta(i) = sum(freq_mutation_in_eta(muts_in_bin)>cutoff)./sum(muts_in_bin);
end
plot(bins+bins(2)/2,(maxyp+2)*prob_in_eta,'k-', 'LineWidth', 2,'MarkerFaceColor', rgb('Black'),'Marker','o')
text(ones(size(0:.2:1)),(maxyp+2)*(0:.2:1),num2str([0:.2:1]'))
plot([.98*ones(1,length(0:.2:1)); ones(1,length(0:.2:1))], (maxyp+2)*[0:.2:1; 0:.2:1], 'k-')



subplot(2,1,2);
hold on;
%calculate
bins=1:max(number_specimens_in);
num_eta=zeros(size(bins));
num_muts=zeros(size(bins));
for i=1:max(number_specimens_in)
    muts_in_bin= number_specimens_in==i;
    num_eta(i) = sum(freq_mutation_in_eta(muts_in_bin)>cutoff);
    num_muts(i) = sum(muts_in_bin);
end
%plot
for i=2:numel(bins)
    patch([bins(i) bins(i) bins(i)+1 bins(i)+1], [0 num_muts(i) num_muts(i) 0], rgb('White'));
    patch([bins(i) bins(i) bins(i)+1 bins(i)+1], [0 num_eta(i) num_eta(i) 0], rgb('Gray'));
end
%make split axis
maxy=max(num_muts(2:end));
maxyp=maxy+30;
ylim([0 maxyp+2])
patch([bins(1) bins(1) bins(2) bins(2)], [0 maxyp maxyp 0], rgb('White'));
patch([bins(1) bins(1) bins(2) bins(2)], [0 num_eta(1) num_eta(1) 0], rgb('Gray'));
set(gca,'Ytick',[0:15:45 maxyp])
set(gca,'Yticklabel',num2str([0:15:45 num_muts(1)]'))
%show probability with split bins for better visualization
prob_in_eta=zeros(numel(bins),1);
num_muts=zeros(numel(bins));
for i=1:numel(bins)
    muts_in_bin=  number_specimens_in > bins(i);
    prob_in_eta(i) = sum(freq_mutation_in_eta(muts_in_bin)>cutoff)./sum(muts_in_bin);
end
plot(bins+bins(1)/2,(maxyp+2)*prob_in_eta,'k-', 'LineWidth', 2,'MarkerFaceColor', rgb('Black'),'Marker','o')
plot([9.8*ones(1,length(0:.2:1)); 10*ones(1,length(0:.2:1))], (maxyp+2)*[0:.2:1; 0:.2:1], 'k-')
text(10*ones(size(0:.2:1)),(maxyp+2)*(0:.2:1),num2str([0:.2:1]'))

xlim([1 10])
set(gca,'Xtick',1.5:1:9.5)
set(gca,'Xticklabel',num2str([1:1:9]'))

%% compare max distance between genotypes, only considering endotreacheal aspirate

d_eta=max_distance_between_genotypes_in_eta(~isnan(max_distance_between_genotypes_in_eta) & max_distance_between_genotypes<100);
d_all=max_distance_between_genotypes(~isnan(max_distance_between_genotypes_in_eta) & max_distance_between_genotypes < 100);
%only including simple infections. mixed strain infection distance isn't
%done properly for aspirate anyway.

% %add some jitter for visibility
jitter=0.5-1*rand(size(d_all));
d_all=d_all+jitter;
d_eta=d_eta+jitter;


epi_threshold = 12;

figure;
clf; hold on;
plot(d_eta,d_all,'ko','MarkerFaceColor',rgb('Gray'),'Markersize', 4)
%plot(0:30,epi_threshold*ones(31,1),'k-')
%plot(epi_threshold*ones(31,1),0:30,'k-')
plot(0:20,0:20,'k:')
set(gca,'Ytick',0:5:20)
set(gca,'Xtick',0:5:20)
ylim([0 16])
xlim([0 16])
xlabel('Max distance between genotypes in eta')
ylabel('Max distance between genotypes in subject')

%% Count transmission events within lung and to spleen and liver, by lung site
%do not include endeotracheal aspirate in either analysis%(doesn't change results significantly)

numtrials=1000;

numtranmissions_to_lung_sites_by_site=nan(numel(subjects),eta_site_number,numtrials);
numtranmissions_lung_to_liver_by_site=nan(numel(subjects),eta_site_number,numtrials);
numtranmissions_lung_to_spleen_by_site=nan(numel(subjects),eta_site_number,numtrials);

for k=1:numel(subjects)
    
    for t=1:numtrials
        
        if ~isempty(genotypematrices{k}) & sum(sum(genotypes_freq_sites_matrices{k}(:,eta_site_number:end),2)>0)>0
            
            %at least 7 samples from each site considered, do not include endotracheal aspirate
            [isolates_per_site, pt_sites] = find_count_duplicates(compartment_lists{k});
            lungsites= pt_sites(pt_sites'<= eta_site_number & isolates_per_site>=min_isolates_per_site);
            
            if numel(lungsites) > 1
                for i=1:numel(lungsites);
                    
                    %get and normalize properly
                    [subsampled_sitetimes, subsampled_sitefreqs] = subsample_from_each_site(compartment_lists{k},observed_matrices_placed_with_multiple_lineages{k},genotypematrices{k},...
                        numel(sites), min_isolates_per_site, snp_threshold_for_assigning_to_genotypes, max_remaining_room_to_add_in_ancestor,max_unexplained_snps_to_add_in_ancestor);
                       
                    %% first do lung analysis
                    
                    %pick a lungsite to remove
                    lungsites_this_interation=lungsites(lungsites~=lungsites(i));
                    sitefreqs=subsampled_sitefreqs;
                    sitefreqs=[sum(sitefreqs(:,lungsites_this_interation),2) sitefreqs(:,lungsites(i))];
                    sitefreqs=sitefreqs./repmat(sum(sitefreqs),size(sitefreqs,1),1);
                    
                    sitetimes=subsampled_sitetimes;
                    sitetimes=[sum(sitetimes(:,lungsites_this_interation),2) sitetimes(:,lungsites(i))];
                    
                    %remove low frequency genotypes from analysis
                    sitefreqs(isnan(sitefreqs))=0;
                    x=sum(sitefreqs(1:end-1,:),2) < min_freq_within_site;
                    sitefreqs(x,:)=[];
                    sitetimes(x,:)=[];
                    genotypes=genotypematrices{k};
                    genotypes(x,:)=[];
                    
                    if ~isempty(genotypes)
                        plotgenotypes=[genotypes; zeros(1,size(genotypes,2))];
                        sitesk=find(sum(sitefreqs)>0);
                        numtranmissions_to_lung_sites_by_site(k,i,t) = count_transmission_on_trees(plotgenotypes, sitetimes, 2, ismember(k,subjects_with_multiple_strains), min_observations_for_transmission_analysis);
                    else
                        numtranmissions_to_lung_sites_by_site(k,i,t) = 1;
                    end
                    
                    
                    %% then go liver/spleen analysis
                    
                    sitefreqs=subsampled_sitefreqs;
                    lungfreqs=sum(sitefreqs(:,lungsites_this_interation),2);
                    sitefreqs(:,1)=lungfreqs;
                    sitefreqs(:,2:eta_site_number)=0;
                    sitefreqs=sitefreqs./repmat(sum(sitefreqs),size(sitefreqs,1),1);
                    
                    sitetimes=subsampled_sitetimes;
                    lungtimes=sum(sitetimes(:,lungsites_this_interation),2);
                    sitetimes(:,1)=lungtimes;
                    sitetimes(:,eta_site_number)=0;
                    
                    %remove low frequency strains from analysis
                    sitefreqs(isnan(sitefreqs))=0;
                    x=sum(sitefreqs(1:end-1,:),2) <min_freq_within_site;
                    sitefreqs(x,:)=[];
                    sitetimes(x,:)=[];
                    
                    genotypes=genotypematrices{k};
                    genotypes(x,:)=[];
                    
                    if ~isempty(genotypes)
                        plotgenotypes=[genotypes; zeros(1,size(genotypes,2))];
                        sitesk=find(sum(sitefreqs)>0);
                        if ismember(liver_site_number,sitesk)
                            numtranmissions_lung_to_liver_by_site(k,i,t) = count_transmission_on_trees(plotgenotypes, sitetimes, liver_site_number, ismember(k,subjects_with_multiple_strains), min_observations_for_transmission_analysis);
                        end
                        if ismember(spleen_site_number,sitesk)
                            numtranmissions_lung_to_spleen_by_site(k,i,t) = count_transmission_on_trees(plotgenotypes, sitetimes, spleen_site_number, ismember(k,subjects_with_multiple_strains), min_observations_for_transmission_analysis);
                        end
                    else
                        if ismember(liver_site_number,sitesk)
                            numtranmissions_lung_to_liver_by_site(k,i,t) = 1;
                        end
                        if ismember(spleen_site_number,sitesk)
                            numtranmissions_lung_to_spleen_by_site(k,i,t) = 1;
                        end
                    end
                end
            end
        end
    end
end





%% compare intralung transmission to interorgan transmission -- Figure 4D


fprintf(fid_F4,'\n\n#Figure 4D\n#Also see Supplementary Figure 4\n');
fprintf(fid_F4,'#Subject, Mean number of transmissions to lung site,Mean number of transmission to lung site (each site - semicolon delimited),Mean number of transmissions to liver,Mean number of transmission to liver (each removed site - semicolon delimited),Mean number of transmissions to spleen,Mean number of transmission to spleen (each removed site - semicolon delimited)\n\n');


figure(28); clf;
hold on;
plot(0:5, 0:5,'-','Color',rgb('Black'))

for k=1:numel(subjects)
    if  sum(~isnan(numtranmissions_to_lung_sites_by_site(k,:,1))) > 0 & ...
            (sum(~isnan(numtranmissions_lung_to_spleen_by_site(k,:,1))) |  sum(~isnan(numtranmissions_lung_to_liver_by_site(k,:,1))))
        
        lungt_mean_each_site=mean(squeeze(numtranmissions_to_lung_sites_by_site(k,:,:)),2);
        lungt_mean_each_site=lungt_mean_each_site(~isnan(lungt_mean_each_site));
        
        spleent_mean_each_site=mean(squeeze(numtranmissions_lung_to_spleen_by_site(k,:,:)),2);
        spleent_mean_each_site=spleent_mean_each_site(~isnan(spleent_mean_each_site));
        
        liver_mean_each_site=mean(squeeze(numtranmissions_lung_to_liver_by_site(k,:,:)),2);
        liver_mean_each_site=liver_mean_each_site(~isnan(liver_mean_each_site));
        
        SEM_lung=std(lungt_mean_each_site)/sqrt(numel(lungt_mean_each_site));
        SEM_liver=std(liver_mean_each_site)/sqrt(numel(liver_mean_each_site));
        SEM_spleen=std(spleent_mean_each_site)/sqrt(numel(spleent_mean_each_site));
        
        if ~isempty(liver_mean_each_site) | ~isempty(spleent_mean_each_site)
            fprintf(fid_F4,[num2str(mean(lungt_mean_each_site)) ',']);
            for i=1:numel(lungt_mean_each_site)
                fprintf(fid_F4,[num2str(lungt_mean_each_site(i)) ';']);
            end
            fprintf(fid_F4,',');
        else
            fprintf(fid_F4,',,');
        end
        
        if ~isempty(liver_mean_each_site)
            jitter=.15-.3*rand(1,1);
            plot([mean(lungt_mean_each_site)+jitter-SEM_lung mean(lungt_mean_each_site)+jitter+SEM_lung], [mean(liver_mean_each_site)+jitter, mean(liver_mean_each_site)+jitter], '-', 'Color', rgb('Gray'),'LineWidth', 1);
            plot([mean(lungt_mean_each_site)+jitter, mean(lungt_mean_each_site)+jitter], [mean(liver_mean_each_site)+jitter-SEM_liver mean(liver_mean_each_site)+jitter+SEM_liver], '-', 'Color', rgb('Gray'),'LineWidth', 1);
            plot(mean(lungt_mean_each_site)+jitter, mean(liver_mean_each_site)+jitter, 'o', 'MarkerFaceColor',  rgb('DeepPink'), 'MarkerEdgeColor', 'none','MarkerSize', 5);
            
            fprintf(fid_F4,[num2str(mean(liver_mean_each_site)) ',']);
            for i=1:numel(liver_mean_each_site)
                fprintf(fid_F4,[num2str(liver_mean_each_site(i)) ';']);
            end
            fprintf(fid_F4,',');
        else
            fprintf(fid_F4,',,');
        end
        
        
        if ~isempty(spleent_mean_each_site)
            jitter=.15-.3*rand(1,1);
            plot([mean(lungt_mean_each_site)+jitter-SEM_lung mean(lungt_mean_each_site)+jitter+SEM_lung], [mean(spleent_mean_each_site)+jitter, mean(spleent_mean_each_site)+jitter], '-', 'Color', rgb('Gray'),'LineWidth', 1);
            plot([mean(lungt_mean_each_site)+jitter, mean(lungt_mean_each_site)+jitter], [mean(spleent_mean_each_site)+jitter-SEM_spleen mean(spleent_mean_each_site)+jitter+SEM_spleen], '-', 'Color', rgb('Gray'),'LineWidth', 1);
            plot(mean(lungt_mean_each_site)+jitter, mean(spleent_mean_each_site)+jitter, 'o', 'MarkerFaceColor',  rgb('Purple'), 'MarkerEdgeColor', 'none', 'MarkerSize', 5);
            
            fprintf(fid_F4,[num2str(mean(spleent_mean_each_site)) ',']);
            for i=1:numel(spleent_mean_each_site)
                fprintf(fid_F4,[num2str(spleent_mean_each_site(i)) ';']);
            end
            fprintf(fid_F4,',');
        else
            fprintf(fid_F4,',,');
        end
        fprintf(fid_F4,'\n');
        
    end
end
xlim([.5 3.5])
ylim([.5 3.5])
set(gca,'Xtick',0:1:4)
set(gca,'Ytick',0:1:4)

xlabel('Number of transmission events to lung specimen')
ylabel('Number of transmission events to liver or spleen')

%%  Make phylogeny figures for transmission -- Figure 4C

fprintf(fid_F4,'\n\n#Figure 4C\n#Also see Supplementary Table 4 for more information on each mutated position and distribution across samples\n');

to_plot_for_transmission = [30 34];

for k=to_plot_for_transmission
    
    fprintf(fid_F4,['\n#Subject P' num2str(k) '\n']);
    fprintf(fid_F4,'#Genomic positions mutated in genotype(semicolon seperated), Number of times observed in lung, Number of times observed in liver, Number of times observed in spleen\n');
    
    genotypes=genotypematrices{k};
    de_novo_muts_pt=de_novo_muts_placed([de_novo_muts_placed(:).subject]==k);
    
    if ~isempty(genotypes) & sum(sum(genotypes_freq_sites_matrices{k}(:,(eta_site_number+1):end),2)>0)>0
        
        %compress and normalize lung
        sitefreqs=genotypes_freq_sites_matrices{k};
        sitetimes=genotype_sites_times_matrices{k};
        
        lungfreqs=sum(sitefreqs(:,1:(eta_site_number)),2);
        sitefreqs(:,1)=lungfreqs;
        sitefreqs(:,2:eta_site_number)=0;
        sitefreqs=sitefreqs./repmat(sum(sitefreqs),size(sitefreqs,1),1);
        sitefreqs(isnan(sitefreqs))=0;
        
        lungtimes=sum(sitetimes(:,1:(eta_site_number)),2);
        sitetimes(:,1)=lungtimes;
        sitetimes(:,2:eta_site_number)=0;
        
        %remove low frequency genotypes from  plotting
        genotypestokeep=sum(sitefreqs(1:end-1,:)>= min_freq_within_site,2) >0 & sum(sitetimes(1:end-1,:),2) >= min_observations_for_transmission_analysis;
        genotypestokeep=[genotypestokeep; true];
        plotsitefreqs = sitefreqs(genotypestokeep,:);
        plotsitetimes = sitetimes(genotypestokeep,:);
        plotgenotypes=[genotypes(genotypestokeep(1:end-1),:); zeros(1,size(genotypes,2))]; %keep ancestor
        
        if ~isempty(genotypes)
            
            %reorder genotypes by number of mutations
            forsorting=sum(plotgenotypes,2)+400*plotgenotypes(:,end);
            [~,sortedorder]=sort(forsorting,'ascend');
            neworder=sortedorder;
            
            if k==30
                neworder=sortedorder([1 2 4 5 6 7 3]);
            elseif k==34
                neworder=sortedorder([1 2 4 5 6 7 3]);
            end
            
            %implement new order
            plotsitetimes=plotsitetimes(neworder,:);
            plotgenotypes=plotgenotypes(neworder,:);
            plotsitefreqs=plotsitefreqs(neworder,:);
            
            %tree
            nts_for_tree=char(plotgenotypes); nts_for_tree(plotgenotypes==0)='A'; nts_for_tree(plotgenotypes==1)='T';
            TreeSampleNames2={};
            TreeColors={};
            for i=1:(size(plotgenotypes,1)+1)
                TreeSampleNames2{end+1}=['GT' num2str(i)];
                TreeColors{end+1}='16777216'; % all same for now
            end
            
            cd([masterdir '/trees/temptreefolder'])
            treefilename=generate_parsimony_tree(nts_for_tree', TreeSampleNames2, ['P' num2str(k) ]);
            color_tree(treefilename, [masterdir '/trees/genotypetrees_organs/P' num2str(k) '.tree'], TreeSampleNames2, TreeColors)
            cd(masterdir)
            
            %make a set of bars for each location
            h=figure; hold on;
            axis off; title(['P' num2str(k)]);
            top=18;
            
            %draw a line for each genotype
            for j=1:size(plotsitefreqs,1);
                plot([0 5], [top-j top-j], 'k:');
            end
            
            sitesk=find(sum(plotsitefreqs)>0);
            for i=1:numel(sitesk)
                genotypesfound=find(plotsitefreqs(:,sitesk(i)));
                freqsincompart=plotsitefreqs(genotypesfound,sitesk(i));
                plotlabels={};
                for j=1:numel(genotypesfound)
                    patch([i i i+.8 i+.8],[top-genotypesfound(j)-.35 top-genotypesfound(j)+.35 top-genotypesfound(j)+.35 top-genotypesfound(j)-.35],...
                        'r','FaceAlpha', 1, 'EdgeColor','none','FaceColor', rgb(site_colors{sitesk(i)}));
                end
            end
            xlim([0 5])
            ylim([0 top+1])
            
            
            %Make CSV file
            for i=1:size(plotsitefreqs,1)
                muts_genotype=find(plotgenotypes(i,:));
                if isempty(muts_genotype)
                    fprintf(fid_F4,'Anc');
                else
                    for j=1:numel(muts_genotype);
                        fprintf(fid_F4,[num2str(de_novo_muts_pt(muts_genotype(j)).pos) ';']);
                    end
                end
                fprintf(fid_F4,',');
                for j=1:numel(sitesk)
                    fprintf(fid_F4,[num2str(plotsitetimes(i,j)) ',']);
                end
                fprintf(fid_F4,'\n');
            end
        end
    end
end