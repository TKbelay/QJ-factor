#run as julia beta_new.jl (#for running on linux command line)
 
using DelimitedFiles, SparseArrays, LinearAlgebra

genot=trunc.(Int64,readdlm("genot_id.ped"))[:];#Read a list of genotyped animal ids
ngenot=trunc.(Int64,readdlm("ngenot_id.ped"))[:]; #Read a list on ungenotyped animal ids

Q_geno=readdlm("ggrpcont_new_genosort_bothid.dat"); #Read genetic group contribution for genotyped animals sorted based renum id
Q2=Q_geno[:,3:end];# the first 2 coulmns might be Id and inbreeding coeffcient

ped=trunc.(Int64,readdlm("renumped_nobd.ped"));# read pedigree file (renumbered)
animals=length(ped[:,1]);

#construct inverse pedigree based relationship matrix following Handersons rule (inbreeding coefficient is not included in the A inv)
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

Ainv=IT'*sparse(Diagonal(ID))*IT;

IA11=Ainv[ngenot,ngenot];# extract block of Ainv for non-genotyped animals
IA12=Ainv[ngenot,genot]; #extract block of Ainv for non-genotyped and genotyped animals

J=-IA11\(IA12*ones(ng)); #calculate J factor, and ng is number of genotyped animals
Qplus=-IA11\(IA12*Q2); #calculate genetic group contributions for non-genotyped animals (Qplus), but orginal genetic contribution (Q2) is used for genotyped animals

writedlm("J_large_valp0.dat",J,"\t")

writedlm("Qplus_matrix.dat",Qplus,"\t")
