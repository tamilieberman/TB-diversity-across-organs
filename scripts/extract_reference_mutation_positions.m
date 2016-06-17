function ref_ID  = extract_reference_mutation_positions(ref_folder, positions)
    
    fastafile = [ref_folder '/genome.fasta']; 
    
   % fprintf('\nGetting outgroup sequences from reference genome\n\t%s\n', fastafile);
    fr = fastaread(fastafile) ;
    
    ref_ID=zeros(size(positions,1),1);
    for i=1:size(positions,1)
        ref_ID(i) = fr(positions(i,1)).Sequence(positions(i,2)); 
    end
    
end