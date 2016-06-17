

%% This part of the pipeline is for identifying de novo mutations within a set of samples


%% display option

show_clickable_tables=0; % set to 1 to make interactive tables for investigating raw data at each variant position


%% define list of subjects

subjects=1:44;
subjects_with_multiple_strains=[11    15    28    37];


%% Parameters  -- IMPORTANT TO ADJUST !!!!
% Much of this will need to vary according to your reference genome,
% coverage, and particular samples

min_average_coverage_to_include_sample = 5;


% for finding fixed mutations between samples
max_fraction_ambigious_samples = .9;
max_mean_coverage_position = 5;
min_qual_for_call = 40;
min_maf_for_call = .7;
min_cov_for_call = 6;
FQ_cutoff=85; %min -FQ that samples supporting both alleles must have


% for finding polymorphic mutations
STRICT_PARAMETERS = struct( 'minorfreqthreshold',           .10, ...
    'maxreads_perstrand_percentile', 99, ...
    'minreads_perstrand',            10, ...
    'minreads_perstrand_per_allele', 3,...
    'min_bq',                        29,...
    'min_mq',                        39, ...
    'min_td',                        25, ... % avg tail dist of each allele
    'max_td',                        75, ...
    'max_sbp',                       3,...  % p-val fisher's exact test of strand bias
    'max_percent_indels',           .20, ...
    'min_control_MAF',              .50, ...
    'max_mqp',                      200, ...% t-test whether map qual dist diff b/t two alleles
    'max_bqp',                      200, ...% t-test whether qual diff between two alleles
    'max_tdp',                      200, ...% t-test whether tail dist diff b/t two alleles
    'max_percent_ends',               1 ...
    );


%for removing positions prone to false-positive polymorphisms
maf_cutoff_for_being_pure = .97;
min_proportion_impure_samples_called_diverse = .15;
min_impure_samples_for_proportion_criteria = 3; %this avoids calling the below criteria for samples with only 10 samples, for which the distribution of mutations across samples is less indicative
max_proprtion_samples_impure = .9;
min_proportion_of_samples_to_be_considered_suspiciously_many = .33; %this goes with next parameter
min_frequency_one_sample_should_be_at_if_many_samples_impure = .6; %this goes with previous parameter

% for lineage-determining positions
cutoff_for_PCA_analysis=.011;

% for removing SNPS due to recombination events
recombination_block_size=1000; %neighborhood in which nearby co-varying SNPs might be due to recominbation event
covariance_cutoff=.5;

% how far upstream of the nearest gene to annotate something a promoter mutation
promotersize=150;

%% Enviornment set up -- probably won't need to change

masterdir=char(pwd);

REFGENOMEFOLDER=[masterdir '/MTB_anc'];
SCRIPTSDIRECTORY = [masterdir '/scripts'];
path(SCRIPTSDIRECTORY,path);

%% Define genomic positions for placing mutations in global phylogeny

% These poisitions, defined in Supplementary Table 3 of
% Coll et al, Nature Communications (2014), Define Global Lineages 1-6
% http://www.nature.com/ncomms/2014/140901/ncomms5812/extref/ncomms5812-s1.pdf

global_lineage_positions=[615938 , 497491 , 3273107 , 931123, 1799921, 1816587];



%% Initialize

NTs='ATCG';

num_samples_removed_based_on_coverage=zeros(size(subjects));


%% Detect mutations within each subject


for k=1:numel(subjects)
    
    
    cd([masterdir '/subject_folders/P' num2str(k)])
    load('candidate_mutation_table')
    
    
    %% Remove undesired samples based on name and/or coverage
    coverage=squeeze(sum(counts(1:8,:,:)));

    Nsample=numel(SampleNames);
    
    Quals = -1*Quals; %use -Quals because this way higher numbers are more confident
    
    
    %% Read in genome information
    
    [ChrStarts, GenomeLength, ~, ScafNames]= genomestats(REFGENOMEFOLDER);
    
    refnt = extract_reference_mutation_positions(REFGENOMEFOLDER, p2chrpos(p,ChrStarts));
    
    
    % Is be more complicated if refnt ~= ancnt
    ancnt=refnt;
    [~,ancnti]=ismember(refnt,NTs);
    ancnti_m=repmat(ancnti,1,Nsample);
    
    annotations = annotate_mutations_gb(p2chrpos(p,ChrStarts),REFGENOMEFOLDER) ;
    
    %% Make some basic structures for finding mutations
    
    [maf, maNT, minorNT, minorAF] = div_major_allele_freq(counts);
    
    [~,refnti]=ismember(refnt,NTs);
    
    mutantAF=zeros(size(maNT));
    mutantAF(maNT~=ancnti_m)=maf(maNT~=ancnti_m);
    mutantAF(minorNT~=ancnti_m)=mutantAF(minorNT~=ancnti_m)+minorAF(minorNT~=ancnti_m); %this construction allows for positions with two different mutations
    
    
    %% Find positions with polymorphic mutations
    
    diversemutation=div_test_thresholds(counts,STRICT_PARAMETERS, coveragethresholds);
    
    %% Find positions with fixed mutations
    
    % Find the nucleotide identity at each position.
    Calls=maNT;
    Calls(Quals < min_qual_for_call | maf< min_maf_for_call | coverage < min_cov_for_call)=0; %Quals < min_qual_for_call |
    Calls(sum(Calls<1,2)>=(Nsample*max_fraction_ambigious_samples) | mean(coverage,2)<max_mean_coverage_position,:)=0;
    
    Calls(diversemutation>0)=-1;  %-1 records that a diverse mutation was called
    
    [MutQual, MutQualIsolates] = ana_mutation_quality(maNT,Quals) ;  %assumes quals already inverted
    fixedmutation=((maNT~=repmat(ancnti,1,Nsample)) & maNT>0 & repmat(MutQual,1,Nsample)>=FQ_cutoff);
    
    
    
    
    
    %% remove adjacent covarying snps which might arise from PPE variation or recombination
    % this approach only makese sense for subjects without multiple
    % lineages, where the expected density of SNPs is low
    % in subjects with multiple lineages, these mutations will be filtered
    % out later, by overall PCA analysis 
    
    
    candidatepositions=find(sum(diversemutation | fixedmutation,2)>0);
        
  %  covariancematrix=cov(mutantAF');
    
    ptoremove=[];
    if ~ismember(k, subjects_with_multiple_strains)
        for j=1:numel(candidatepositions);
            i=candidatepositions(j);
            %find nearby snps
            region=find(p>p(i)-recombination_block_size,1):find(p<p(i)+recombination_block_size,1,'last');
            if numel(region)>1
                %for each pair in here, if high enough covariance, discard
                covariancematrix = cov(mutantAF(region,:)');
                max_covariance = covariance_cutoff * max(diag(covariancematrix),[],2); 
                [a,b]=find(covariancematrix>repmat(max_covariance,1,length(covariancematrix)) & covariancematrix > 0.001);
                ptoremove=[ptoremove region(a(a~=b)) region(b(a~=b))];
            end
        end
        
        ptoremove=unique(ptoremove);
    end
                 
    
    diversemutation(ptoremove,:)=0;
    fixedmutation(ptoremove,:)=0;
    
    
    
    %%
    
    
    hasmutation_all= fixedmutation | diversemutation;
    
    num_impure_samples = sum(maf<maf_cutoff_for_being_pure,2);
    
    %Remove troublesome genomic positions
    bad_positions_diverse_strict= (num_impure_samples >= floor(max_proprtion_samples_impure*Nsample)) | ...
        sum(maf<maf_cutoff_for_being_pure,2) >= min_impure_samples_for_proportion_criteria &...
        (sum(hasmutation_all>0,2)./num_impure_samples <= min_proportion_impure_samples_called_diverse |...
        (num_impure_samples >= min_proportion_of_samples_to_be_considered_suspiciously_many*Nsample) & sum(mutantAF>min_frequency_one_sample_should_be_at_if_many_samples_impure,2)<1);
    %too many impure samples OR
    %many impure strains, but too few passed filters OR
    %many impure mutations, but none with mutation at high frequency
    
    
    
    %this less strict threshold important for discovering mutations
    %betweeen strains, when calculating dMRCA
    bad_positions_diverse_lessstrict= (num_impure_samples >= floor(max_proprtion_samples_impure*Nsample)) | ...
        num_impure_samples>min_impure_samples_for_proportion_criteria &...
        sum(hasmutation_all>0,2)./num_impure_samples < min_proportion_impure_samples_called_diverse;
    %too many impure samples OR
    %many impure strains, but too few passed filters
    
    
    %% find lineage determining mutations
    
    % find lineage determinents
    
    
    
    
    %use stricter threshold for finding de novo SNPS and illustration in
    %Supp figure.
    %lesser threshold includes false positives that screw up correlation,
    %but these only minimally inflate the calculation of <dMRCA> but leaving them out
    %greatly depresses it for Subject 757
    
    
    
    
    hasmutation_for_lineages=hasmutation_all; hasmutation_for_lineages(bad_positions_diverse_strict,:)=0;
    goodindexfor_lineages=find(sum(hasmutation_for_lineages,2)>0);
    mutantAF_for_lineages=mutantAF(goodindexfor_lineages,:); mutantAF_for_lineages(isnan(mutantAF_for_lineages))=0;
    [coeff, score, latent, tsquared, explained] = pca(mutantAF_for_lineages');
    
    %make a plot about it
    if size(coeff,2)>1 & ismember(k,subjects_with_multiple_strains)
        lineagedeterminents=goodindexfor_lineages(abs(coeff(:,1))>cutoff_for_PCA_analysis);
        figure;clf; subplot(1, 2, 1);hold on;
        title(['P' num2str(k)])
        [bins,patches]=hist(coeff(:,1),-.05:.005:.05);
        plot([-cutoff_for_PCA_analysis -cutoff_for_PCA_analysis], [0 500], ':', 'LineWidth', 1, 'Color', rgb('Black'))
        plot([cutoff_for_PCA_analysis cutoff_for_PCA_analysis], [0 500], ':', 'LineWidth', 1, 'Color', rgb('Black'))
        bar(patches(patches< -cutoff_for_PCA_analysis), bins(patches< -cutoff_for_PCA_analysis), 'FaceColor', rgb('DarkGrey'),'EdgeColor','none')
        bar(patches(patches > cutoff_for_PCA_analysis), bins(patches > cutoff_for_PCA_analysis), 'FaceColor', rgb('LightGrey'),'EdgeColor','none')
        bar(patches(patches > -cutoff_for_PCA_analysis & patches < cutoff_for_PCA_analysis), bins(patches > -cutoff_for_PCA_analysis & patches < cutoff_for_PCA_analysis), 'FaceColor', rgb('Magenta'),'EdgeColor','none')
        ylim([0 max(bins)*1.05])
        xlim([-.06 .06])
        set(gca,'Ytick',0:50:200)
        set(gca,'Xtick',[-.05 0 .05])
        xlabel('Coefficient of mutation in first principal component');
        ylabel('Number of mutations')
        
        
        subplot(1, 2, 2); hold on;
        title(['P' num2str(k)])
        lineage2determinents=goodindexfor_lineages(coeff(:,1)>cutoff_for_PCA_analysis);
        lineage1determinents=goodindexfor_lineages(coeff(:,1)<-cutoff_for_PCA_analysis);
        dennovomutp=goodindexfor_lineages(abs(coeff(:,1))<cutoff_for_PCA_analysis);
        x=mutantAF(lineage2determinents,:); x(isnan(x))=0;
        proportionlineage2=mean(x);
        [bins,patches]=hist(1-proportionlineage2,10);
        bar(patches,bins,'FaceColor', rgb('DarkGrey'),'EdgeColor','none')
        xlabel('Frequency of Strain 1 in sample');
        set(gca,'Ytick',0:20:60)
        ylabel('Number of samples')
        xlim([-.01 1.01])
        
        %%
        if k==15
            %Plot for Supplementary Figure 1c
            %samples=datasample(1:numel(SampleNames),5,'Replace', false);
            samples=[88 53 64 31 26];
             clf; hold on;
            for i=1:numel(samples)
                bins=0:.05:1;
                subplot(numel(samples),1,i)
                x=mutantAF(lineage2determinents,samples(i));
                [bins2,patches2]=hist(x(~isnan(x)),bins);
                x=mutantAF(lineage1determinents,samples(i));
                [bins1,patches1]=hist(x(~isnan(x)),bins);
                x=mutantAF(dennovomutp,samples(i));
                [bins3,patches3]=hist(x(~isnan(x)),bins);
                bar(patches1,[bins1; bins2; bins3]','Stacked','EdgeColor','none')
                colormap([rgb('LightGrey'); rgb('DarkGrey'); rgb('Magenta')])
                xlim([0 1])
            end
        end

        %[88 53 64 31 80]
        
    else
        lineagedeterminents=[];
    end
    
    
    
    
    %% redo entire calculation with less strict tresholds for better <dMRCA> estimation
    if  ismember(k,subjects_with_multiple_strains)
        hasmutation_for_lineages2=hasmutation_all; hasmutation_for_lineages2(bad_positions_diverse_lessstrict,:)=0;
        goodindexfor_lineages2=find(sum(hasmutation_for_lineages2,2)>0);
        x=mutantAF(goodindexfor_lineages2,:); x(isnan(x))=0;
        [coeff2, ~, ~, ~, ~] = pca(x');
        num_muts_seperating_strains=[sum(coeff2(:,1)<-cutoff_for_PCA_analysis) sum(coeff2(:,1)>cutoff_for_PCA_analysis)];
    end
    
    %% remove bad postions
    
    diversemutation(bad_positions_diverse_strict,:)=0;
    fixedmutation(bad_positions_diverse_strict,:)=0;
    %  diversemutation(ismember(p,problematic_positions_on_reference_genome),:)=0; %manually done to remove false positives
    
    hasmutation=fixedmutation | diversemutation;
    
    %% global lineages
    
    %based on Coll et al 2014, http://www.nature.com/ncomms/2014/140901/ncomms5812/extref/ncomms5812-s1.pdf
    
    proportion_global_lineages=zeros(numel(global_lineage_positions),numel(SampleNames));
    [glp_found,glp_position]=ismember(global_lineage_positions,p);
    glp_found=find(glp_found);
    for i=1:numel(glp_found)
        proportion_global_lineages(glp_found(i),:)=mean(mutantAF(glp_position(glp_found(i)),:));
    end
    
    %% remove lineage defining mutations for de novo muts
    
    order=1:numel(SampleNames);
    
    
    goodpos=sum(hasmutation,2)>0;
    has_denovo_mutation=hasmutation;
    
    lineagedeterminents=ismember(1:numel(p),lineagedeterminents);
    if ismember(k, subjects_with_multiple_strains)
        multiplelineages=1;
        lineage1determinents=goodindexfor_lineages(coeff(:,1)<-cutoff_for_PCA_analysis);
        x=mutantAF(lineage1determinents,:); x(isnan(x))=0;
        proportionlineage1=mean(x);
        has_denovo_mutation(lineagedeterminents,:)=0;
        goodpos(lineagedeterminents)=0;
    else
        multiplelineages=0;
        proportionlineage1=ones(size(SampleNames));
        num_muts_seperating_strains=0;
    end
    
    
    
    %%
    
    annotation_full= append_annotations(annotations, ancnti, Calls, counts, has_denovo_mutation, promotersize) ; %% adds information about particular mutations observed, based on
    
    QualSort=0; %set 1 to show mutations with lowest FQ scores up top, 0 to show in order on the genome
    
    if show_clickable_tables==1;
        clickable_snp_table(annotation_full(goodpos), Calls(goodpos,order), counts(:,goodpos,order), SampleNames(order), ScafNames, MutQual(goodpos), QualSort);
    end
    
    %%  save table
    
    de_novo_muts=annotation_full(goodpos);
    averagecoverage=coveragethresholds(50,:);
    observedmatrix=mutantAF(goodpos,:);
    
    if numel(de_novo_muts) > 50
        fprintf(1,['Caution: Sanity check -- detected ' num2str(numel(de_novo_muts)) ' de novo mutations in P' num2str(k)]);
    end
    
    if ~isfield(de_novo_muts,'muts') %make sure field is there
        for i=1:numel(de_novo_muts)
            de_novo_muts(i).muts={};
        end
    end
    
    save('de_novo_muts','de_novo_muts','SampleNames','observedmatrix', 'averagecoverage', 'proportionlineage1', 'num_muts_seperating_strains','proportion_global_lineages')
    fprintf(1,['\nSubject: P' num2str(k)]);
    fprintf(1,['\nDe novo muts: ' num2str(numel(de_novo_muts))]);
    if ismember(k,subjects_with_multiple_strains)
        fprintf(1,['\nNumber of lineage-defining mutations: ' num2str(num_muts_seperating_strains)]);
    end
    fprintf(1,'\nProportion each global lineage:\n');
    disp(mean(proportion_global_lineages,2))
    
end


