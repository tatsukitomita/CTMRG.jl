using Test
import CTMRG
using CTMRG: ising_transfer_element_1D, ising_transfer_matrix_1D, ising_transfer_element_2D, ising_transfer_tensor_2D
using CTMRG: make_inner_tensors, contruct_with_boundary
using CTMRG: make_2D_index_σtags, make_2D_index_stags
using CTMRG: make_diag_Tesnor
using CTMRG: Index, inds

@testset "CTMRG.jlのテスト" begin
    #１次元のイジング模型の分配関数を計算する関数のテスト
    @testset "ising_transfer_element_1D" begin
        β = 1.0
        J = 1.0
        @test ising_transfer_element_1D(0, 0, β, J) ≈ exp(β*J)
        @test ising_transfer_element_1D(0, 1, β, J) ≈ exp(-β*J)
        @test ising_transfer_element_1D(1, 0, β, J) ≈ exp(-β*J)     
        @test ising_transfer_element_1D(1, 1, β, J) ≈ exp(β*J)    
    end
        @testset "ising_transfer_matrix_1Dのテスト" begin
        β = 1.0
        J = 1.0
        s1 = Index(2, "s1,m=1")
        s2 = Index(2, "s2,n=1")
        T = ising_transfer_matrix_1D(β,J,s1,s2)
        @test T[s1=>1, s2=>1] ≈ exp(β*J)
    end
    #2次元イジング模型の分配関数を計算する関数のテスト
    @testset "ising_transfer_element_2Dのテスト" begin
        β = 1.0
        J = 1.0
        @test ising_transfer_element_2D(1, 1, 1, 1, β, J) ≈ exp(β*J*4) #全てのスピンがup
        @test ising_transfer_element_2D(0, 0, 0, 0, β, J) ≈ exp(β*J*4) #全てのスピンがdown
        @test ising_transfer_element_2D(1, 0, 0, 0, β, J) ≈ 1.0
        @test ising_transfer_element_2D(1, 1, 0, 0, β, J) ≈ exp(-β*J*4)
    end
    @testset "ising_transfer_tensor_2Dのテスト" begin
        σ1 = Index(2, "σ1,m=1")
        σ2 = Index(2, "σ2,m=1")
        s1 = Index(2, "s1,n=1")
        s2 = Index(2, "s2,n=1")
        β = 1.0
        J = 1.0
        T = ising_transfer_tensor_2D(β, J, σ1, σ2, s1, s2)
        @test T[σ1=>1, σ2=>1, s1=>1, s2=>1] ≈ exp(β*J*4) #全てのスピンがup
        @test T[σ1=>2, σ2=>2, s1=>2, s2=>2] ≈ exp(β*J*4) #全てのスピンがdown
        @test T[σ1=>1, σ2=>2, s1=>2, s2=>2] ≈ 1.0
        @test T[σ1=>1, σ2=>1, s1=>2, s2=>2] ≈ exp(-β*J*4)
    end
    @testset "make_inner_tensorsのテスト" begin
        σ = CTMRG.make_2D_index_σtags(1)
        s = CTMRG.make_2D_index_stags(1)
        T_list = make_inner_tensors(1, β, J, s, σ)
        for i in 1:size(T_list, 1), j in 1:size(T_list, 2)
            @test T_list[i,j][σ[j,i]=>1, σ[j+1,i]=>1, s[j,i]=>1, s[j,i+1]=>1] ≈ exp(β*J*4) #全てのスピンがup
            @test T_list[i,j][σ[j,i]=>2, σ[j+1,i]=>2, s[j,i]=>2, s[j,i+1]=>2] ≈ exp(β*J*4) #全てのスピンがdown
            @test T_list[i,j][σ[j,i]=>1, σ[j+1,i]=>2, s[j,i]=>2, s[j,i+1]=>2] ≈ 1.0
            @test T_list[i,j][σ[j,i]=>1, σ[j+1,i]=>1, s[j,i]=>2, s[j,i+1]=>2] ≈ exp(-β*J*4)
        end
    end
    @testset "contract_with_boundaryのテスト" begin
        L = 4
        β = 1.0
        J = 1.0

        T_list, T1_list = contruct_with_boundary(L, β, J)

        @test size(T_list) == size(T1_list) == (2L, 2L)

        for i in 1:2L, j in 1:2L
            inds_orig = Set(inds(T_list[i,j]))
            inds_new  = Set(inds(T1_list[i,j]))

            is_boundary = (i == 1 || i == 2L || j == 1 || j == 2L)

            if is_boundary
                @test inds_orig != inds_new  # 境界のテンソルは ind が変化しているはず
            else
                @test inds_orig == inds_new  # 内部のテンソルは変化していないはず
            end
        end
    end
    #diag.jlのテスト
    @testset "make_diag_Tensorのテスト" begin
        new_1 = Index(2, "new = 1")
        new_2 = Index(2, "new = 2")
        new1_list =[new_1]
        new2_list = [new_2]
        C = [2.0 1.0; 1.0 2.0]
        evals, evecs = eigen(C)
        ievecs = inv(evecs)
        Λ, A, iA = make_diag_Tesnor(evals, evecs, ievecs)
    end

end