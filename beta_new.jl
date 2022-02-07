#run as julia beta_new.jl (#for julia version 1.4.2)
 
using DelimitedFiles
genot=trunc.(Int64,readdlm("genot_id.ped"))[:];#Read a list of genotyped animal id for 30729
ngenot=trunc.(Int64,readdlm("ngenot_id.ped"))[:]; #Read a list on ungenotyped animal id

Q_geno=readdlm("ggrpcont_new_genosort_bothid.dat"); #Read genetic group contribution for genotyped animals sorted based renum id
Q=Q_geno[:,3:end];#307296x115

ped=trunc.(Int64,readdlm("renumped_nobd.ped"));
#ped=trunc.(Int64,readdlm("renumped_valp0_nobd.ped"));
animals=length(ped[:,1]);

using SparseArrays
using LinearAlgebra
IT=sparse(1.0*I,animals,animals);#calculate inverse of the lower triangular matrix (T) in A-matrix: (inv(A)=(inv(T))'*inv(D)*inv(T)
ID=2*ones(animals);#calculate inverse of diagonal matrix(D), which is reciprocal of diagonal elements of A-matrix

for i=1:animals
    far=ped[i,2];
    mor=ped[i,3];
    if (far*mor)==0
        ID[i]=4/3;
    end
    if (far+mor)==0
       ID[i]=1;
    end
#elements in inverse D are finshed;
   if far!=0
      IT[i,far]=-0.5;
   end
   if mor!=0
      IT[i,mor]=-0.5;
   end
#elements in inverse T are finshed;
end

using LinearAlgebra

Ainv=IT'*sparse(Diagonal(ID))*IT;

IA11=Ainv[ngenot,ngenot];
IA12=Ainv[ngenot,genot];

J=-IA11\(IA12*ones(30729));
Qplus=-IA11\(IA12*Q);

writedlm("J_large_valp0.dat",J,"\t")

writedlm("Qplus_matrix.dat",Qplus,"\t")

#writedlm("Qplusnew_matrix_valp0.dat",J_Q,"\t")
