using Test, BioSequences


@testset "Encoding" begin
    seq = LongDNA{4}("ACGTA")
    @test wgregseq.utils.onehot_encoder(seq) == [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 0 0 0]
    seq = LongDNA{4}("ACGTN")
    @test_throws ArgumentError wgregseq.utils.onehot_encoder(seq)
end

@testset "Energy Matrix" begin
    seq1 = LongDNA{4}("ACGTAA")
    seq2 = LongDNA{4}("AAAAAA")
    emat = zeros(6, 4)
    emat[:, 1] .= 1
    
    @test wgregseq.utils.eval_emat(seq1, emat) == 3
    @test wgregseq.utils.eval_emat(seq2, emat) == 6
    
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, zeros(5, 4))
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, zeros(6, 3))
    
    emat = zeros(4, 4)
    emat[:, 1] .= 1
    @test wgregseq.utils.eval_emat(seq2, emat, start_ind_seq=1, stop_ind_seq=4) == 4
    
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_seq=1)
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_seq=-1, stop_ind_seq=4)
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_seq=-1, stop_ind_seq=4)
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_seq=3, stop_ind_seq=4)
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_seq=3, stop_ind_seq=4)
    
    emat = zeros(6, 4)
    emat[:, 1] .= 1
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_emat=-1, stop_ind_emat=4)
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_emat=-1, stop_ind_emat=4)
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_emat=3, stop_ind_emat=4)
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_emat=3, stop_ind_emat=4)
    
    @test wgregseq.utils.eval_emat(seq2, emat, start_ind_seq=3, stop_ind_seq=4, start_ind_emat=3, stop_ind_emat=4) == 2
    @test_throws ArgumentError wgregseq.utils.eval_emat(seq2, emat, start_ind_seq=3, stop_ind_seq=4, start_ind_emat=2, stop_ind_emat=4)
    
end
