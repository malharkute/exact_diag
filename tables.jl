# Reference tables used for multiplet codes, can be found in TPD's Mutliplet Tutorial
#These tables often use 0-index when indexing with orbital angular momentum (l)
#These tables often use -l:l indexes when indexing with magnetic quantum number (m_l)

using OffsetArrays

#Examples for using OffsetArrays:
#array_name = OffsetVector([element1, element2, element3, element3], start_index:end_index)
#array_name = OffsetArray{Float64}(undef, -1:1, -7:7)

function allowed_k_table()
    # array of allowed k values (gaunt coefficient k) for different transitions combinations
    # indexing is 0:3 (for the orbital angular momentum Q#s)
    
    k_table = OffsetArray{Vector{Int16}}(undef,0:3,0:3);
    
    # s row
    k_table[0,0] = [0]
    k_table[0,1] = [1]
    k_table[0,2] = [0,2]
    k_table[0,3] = [1,3]
    
    # p row
    k_table[1,0] = [1]
    k_table[1,1] = [0,2]
    k_table[1,2] = [1,3]
    k_table[1,3] = [0,2,4]
    
    # d row
    k_table[2,0] = [0,2]
    k_table[2,1] = [1,3]
    k_table[2,2] = [0,2,4]
    k_table[2,3] = [1,3,5]
    
    # f row
    k_table[3,0] = [1,3]
    k_table[3,1] = [0,2,4]
    k_table[3,2] = [1,3,5]
    k_table[3,3] = [0,2,4,6]
    
    k_table
end

function gaunt_coeff_table()
    # This function will generate the Gaunt coefficient table
    # orbital angular momentum indexing is 0:3 for l
    # second indexing for magnetic quantum number is -l:l for m_l
    # (negative values are for change in sign between m1/m2- used for first ml indexing only
    
    #Example:
    #table[l1][12][m1][m2] = [array of gaunt coefficients] <- like table 9/10 in TPD's Multiplet Tutorial
    
    #NOTE: l1 < l2 only!!! Like in TPD's table
    
    gaunt_table = OffsetArray{OffsetArray{OffsetArray{OffsetArray{Vector{Float64}}}}}(undef, 0:3)
    
    # =============== l1 = s set up ==================
    gaunt_table[0] = OffsetArray{OffsetArray{OffsetArray{Vector{Float64}}}}(undef, 0:3)
    
    # s-s transitions
    gaunt_table[0][0] = OffsetArray{OffsetArray{Vector{Float64}}}(undef, 0:0)

    gaunt_table[0][0][0] = OffsetArray{Vector{Float64}}(undef, 0:0)
    gaunt_table[0][0][0][0] = [1,0]
    
    # s-p transitions
    gaunt_table[0][1] = OffsetArray{OffsetArray{Vector{Float64}}}(undef, 0:0)

    gaunt_table[0][1][0] = OffsetArray{Vector{Float64}}(undef, -1:1)
    gaunt_table[0][1][0][1] = [-1,0]
    gaunt_table[0][1][0][0] = [1,0]
    gaunt_table[0][1][0][-1] = [-1,0]
    
    # s-d transitions
    gaunt_table[0][2] = OffsetArray{OffsetArray{Vector{Float64}}}(undef, 0:0)
    
    gaunt_table[0][2][0] = OffsetArray{Vector{Float64}}(undef, -2:2)
    gaunt_table[0][2][0][2] = [0,1]
    gaunt_table[0][2][0][1] = [0,-1]
    gaunt_table[0][2][0][0] = [0,1]
    gaunt_table[0][2][0][-1] = [0,-1]
    gaunt_table[0][2][0][-2] = [0,1]
    
    #s-f transitions
    gaunt_table[0][3] = OffsetArray{OffsetArray{Vector{Float64}}}(undef, 0:0)
    
    gaunt_table[0][3][0] = OffsetArray{Vector{Float64}}(undef, -3:3)
    gaunt_table[0][3][0][3] = [-1]
    gaunt_table[0][3][0][2] = [1]
    gaunt_table[0][3][0][1] = [-1]
    gaunt_table[0][3][0][0] = [1]
    gaunt_table[0][3][0][-1] = [-1]
    gaunt_table[0][3][0][-2] = [1]
    gaunt_table[0][3][0][-3] = [-1]
    
    # =============== l1 = p set up ==================
    gaunt_table[1] = OffsetArray{OffsetArray{OffsetArray{Vector{Float64}}}}(undef, 1:3)
    
    #p-p transitions
    gaunt_table[1][1] = OffsetArray{OffsetArray{Vector{Float64}}}(undef, -1:1)
    
    gaunt_table[1][1][1] = OffsetArray{Vector{Float64}}(undef, -1:1)
    gaunt_table[1][1][1][1] = [1,-1]
    gaunt_table[1][1][1][0] = [0, sqrt(3)]
    gaunt_table[1][1][1][-1] = [0, -1*sqrt(6)]
    
    gaunt_table[1][1][0] = OffsetArray{Vector{Float64}}(undef, 0:0)
    gaunt_table[1][1][0][0] = [1,2]
    
    gaunt_table[1][1][-1] = OffsetArray{Vector{Float64}}(undef, -1:1)
    gaunt_table[1][1][-1][-1] = [1,-1]
    gaunt_table[1][1][-1][0] = [0, sqrt(3)]
    gaunt_table[1][1][-1][1] = [0, -1*sqrt(6)]
    
    #p-d transitions
    gaunt_table[1][2] = OffsetArray{OffsetArray{Vector{Float64}}}(undef, -1:1)
    
    gaunt_table[1][2][1] = OffsetArray{Vector{Float64}}(undef, -2:2)
    gaunt_table[1][2][1][2] = [-1*sqrt(6), sqrt(3)]
    gaunt_table[1][2][1][1] = [sqrt(3), -3]
    gaunt_table[1][2][1][0] = [-1, sqrt(18)]
    gaunt_table[1][2][1][-1] = [0, -1*sqrt(30)]
    gaunt_table[1][2][1][-2] = [0, sqrt(45)]
    
    gaunt_table[1][2][0] = OffsetArray{Vector{Float64}}(undef, -2:2)
    gaunt_table[1][2][0][2] = [0, sqrt(15)]
    gaunt_table[1][2][0][1] = [-1*sqrt(3), -1*sqrt(24)]
    gaunt_table[1][2][0][0] = [2, sqrt(27)]
    gaunt_table[1][2][0][-2] = [0, sqrt(15)]
    gaunt_table[1][2][0][-1] = [-1*sqrt(3), -1*sqrt(24)]
    
    gaunt_table[1][2][-1] = OffsetArray{Vector{Float64}}(undef, -2:2)
    gaunt_table[1][2][-1][-2] = [-1*sqrt(6), sqrt(3)]
    gaunt_table[1][2][-1][-1] = [sqrt(3), -3]
    gaunt_table[1][2][-1][0] = [-1, sqrt(18)]
    gaunt_table[1][2][-1][1] = [0, -1*sqrt(30)]
    gaunt_table[1][2][-1][2] = [0, sqrt(45)]
    
    #p-f transitions
    #COMING SOON - see Condon, E.U & Shortley, G.H. (1953). The Theory of Atomic Spectra. Cambridge.
    
    # =============== l1 = d set up ==================
    gaunt_table[2] = OffsetArray{OffsetArray{OffsetArray{Vector{Float64}}}}(undef, 2:3)
    
    #d-d transitions
    gaunt_table[2][2] = OffsetArray{OffsetArray{Vector{Float64}}}(undef, -2:2)
    
    gaunt_table[2][2][2] = OffsetArray{Vector{Float64}}(undef, -2:2)
    gaunt_table[2][2][2][2] = [1,-2,1]
    gaunt_table[2][2][2][1] = [0,sqrt(6),-1*sqrt(5)]
    gaunt_table[2][2][2][0] = [0,-2,sqrt(15)]
    gaunt_table[2][2][2][-1] = [0,0,-1*sqrt(35)]
    gaunt_table[2][2][2][-2] = [0,0,sqrt(70)]
    
    gaunt_table[2][2][1] = OffsetArray{Vector{Float64}}(undef, -1:1)
    gaunt_table[2][2][1][1] = [1,1,-4]
    gaunt_table[2][2][1][0] = [0,1,sqrt(30)]
    gaunt_table[2][2][1][-1] = [0,-1*sqrt(6),-1*sqrt(40)]
    
    gaunt_table[2][2][0] = OffsetArray{Vector{Float64}}(undef, 0:0)
    gaunt_table[2][2][0][0] = [1,2,6]
    
    gaunt_table[2][2][-1] = OffsetArray{Vector{Float64}}(undef, -1:1)
    gaunt_table[2][2][-1][-1] = [1,1,-4]
    gaunt_table[2][2][-1][0] = [0,1,sqrt(30)]
    gaunt_table[2][2][-1][1] = [0,-1*sqrt(6),-1*sqrt(40)]
    
    gaunt_table[2][2][-2] = OffsetArray{Vector{Float64}}(undef, -2:2)
    gaunt_table[2][2][-2][-2] = [1,-2,1]
    gaunt_table[2][2][-2][-1] = [0,sqrt(6),-1*sqrt(5)]
    gaunt_table[2][2][-2][0] = [0,-2,sqrt(15)]
    gaunt_table[2][2][-2][1] = [0,0,-1*sqrt(35)]
    gaunt_table[2][2][-2][2] = [0,0,sqrt(70)]
    
    #d-f transitions
    gaunt_table[2][3] = OffsetArray{OffsetArray{Vector{Float64}}}(undef, -2:2)
    
    gaunt_table[2][3][2] = OffsetArray{Vector{Float64}}(undef, -3:3)
    gaunt_table[2][3][2][3] = [-1*sqrt(15),sqrt(10),-1]
    gaunt_table[2][3][2][2] = [sqrt(5),-2*sqrt(5),sqrt(5)]
    gaunt_table[2][3][2][1] = [-1,2*sqrt(6),-1*sqrt(15)]
    gaunt_table[2][3][2][0] = [0,-2*sqrt(5),sqrt(35)]
    gaunt_table[2][3][2][-1] = [0,sqrt(10),-1*sqrt(70)]
    gaunt_table[2][3][2][-2] = [0,0,3*sqrt(14)]
    gaunt_table[2][3][2][-3] = [0,0,-1*sqrt(210)]
    
    gaunt_table[2][3][1] = OffsetArray{Vector{Float64}}(undef, -3:3)
    gaunt_table[2][3][1][3] = [0,5,-1*sqrt(7)]
    gaunt_table[2][3][1][2] = [-1*sqrt(10),-1*sqrt(15),2*sqrt(6)]
    gaunt_table[2][3][1][1] = [2*sqrt(2),sqrt(2),-5*sqrt(2)]
    gaunt_table[2][3][1][0] = [-1*sqrt(3),sqrt(2),4*sqrt(5)]
    gaunt_table[2][3][1][-1] = [0,-1*sqrt(15),-1*sqrt(105)]
    gaunt_table[2][3][1][-2] = [0,5,4*sqrt(7)]
    gaunt_table[2][3][1][-3] = [0,0,-2*sqrt(21)]
    
    gaunt_table[2][3][0] = OffsetArray{Vector{Float64}}(undef, -3:3)
    gaunt_table[2][3][0][3] = [0,5,-2*sqrt(7)]
    gaunt_table[2][3][0][2] = [0,0,3*sqrt(7)]
    gaunt_table[2][3][0][1] = [-1*sqrt(6),-3,-3*sqrt(10)]
    gaunt_table[2][3][0][0] = [3,4,10]
    gaunt_table[2][3][0][-3] = [0,5,-2*sqrt(7)]
    gaunt_table[2][3][0][-2] = [0,0,3*sqrt(7)]
    gaunt_table[2][3][0][-1] = [-1*sqrt(6),-3,-3*sqrt(10)]
    
    gaunt_table[2][3][-1] = OffsetArray{Vector{Float64}}(undef, -3:3)
    gaunt_table[2][3][-1][-3] = [0,5,-1*sqrt(7)]
    gaunt_table[2][3][-1][-2] = [-1*sqrt(10),-1*sqrt(15),2*sqrt(6)]
    gaunt_table[2][3][-1][-1] = [2*sqrt(2),sqrt(2),-5*sqrt(2)]
    gaunt_table[2][3][-1][0] = [-1*sqrt(3),sqrt(2),4*sqrt(5)]
    gaunt_table[2][3][-1][1] = [0,-1*sqrt(15),-1*sqrt(105)]
    gaunt_table[2][3][-1][2] = [0,5,4*sqrt(7)]
    gaunt_table[2][3][-1][3] = [0,0,-2*sqrt(21)]
    
    gaunt_table[2][3][-2] = OffsetArray{Vector{Float64}}(undef, -3:3)
    gaunt_table[2][3][-2][-3] = [-1*sqrt(15),sqrt(10),-1]
    gaunt_table[2][3][-2][-2] = [sqrt(5),-2*sqrt(5),sqrt(5)]
    gaunt_table[2][3][-2][-1] = [-1,2*sqrt(6),-1*sqrt(15)]
    gaunt_table[2][3][-2][0] = [0,-2*sqrt(5),sqrt(35)]
    gaunt_table[2][3][-2][1] = [0,sqrt(10),-1*sqrt(70)]
    gaunt_table[2][3][-2][2] = [0,0,3*sqrt(14)]
    gaunt_table[2][3][-2][3] = [0,0,-1*sqrt(210)]
    
    # =============== l1 = f set up ==================
    gaunt_table[3] = OffsetArray{OffsetArray{OffsetArray{Vector{Float64}}}}(undef, 3:3)
    
    #f-f transitions
    gaunt_table[3][3] = OffsetArray{OffsetArray{Vector{Float64}}}(undef, -3:3)
    
    gaunt_table[3][3][3] = OffsetArray{Vector{Float64}}(undef, -3:3)
    gaunt_table[3][3][3][3] = [1,-5,3,-1]
    gaunt_table[3][3][3][2] = [0,5,-1*sqrt(30),sqrt(7)]
    gaunt_table[3][3][3][1] = [0,-1*sqrt(10),sqrt(54),-1*sqrt(28)]
    gaunt_table[3][3][3][0] = [0,0,-1*sqrt(63),sqrt(84)]
    gaunt_table[3][3][3][-1] = [0,0,sqrt(42),-1*sqrt(210)]
    gaunt_table[3][3][3][-2] = [0,0,0,sqrt(462)]
    gaunt_table[3][3][3][-3] = [0,0,0,-1*sqrt(924)]
    
    gaunt_table[3][3][2] = OffsetArray{Vector{Float64}}(undef, -2:2)
    gaunt_table[3][3][2][2] = [1,0,-7,6]
    gaunt_table[3][3][2][1] = [0,sqrt(15),sqrt(32),-1*sqrt(105)]
    gaunt_table[3][3][2][0] = [0,-1*sqrt(20),-1*sqrt(3),4*sqrt(14)]
    gaunt_table[3][3][2][-1] = [0,0,-1*sqrt(14),-1*sqrt(378)]
    gaunt_table[3][3][2][-2] = [0,0,sqrt(70),sqrt(504)]
    
    gaunt_table[3][3][1] = OffsetArray{Vector{Float64}}(undef, -1:1)
    gaunt_table[3][3][1][1] = [1,3,1,-15]
    gaunt_table[3][3][1][0] = [0,sqrt(2),sqrt(15),5*sqrt(14)]
    gaunt_table[3][3][1][-1] = [0,-1*sqrt(24),-1*sqrt(40),-1*sqrt(420)]
    
    gaunt_table[3][3][0] = OffsetArray{Vector{Float64}}(undef, 0:0)
    gaunt_table[3][3][0][0] = [1,4,6,20]
    
    gaunt_table[3][3][-1] = OffsetArray{Vector{Float64}}(undef, -1:1)
    gaunt_table[3][3][-1][-1] = [1,3,1,-15]
    gaunt_table[3][3][-1][0] = [0,sqrt(2),sqrt(15),5*sqrt(14)]
    gaunt_table[3][3][-1][1] = [0,-1*sqrt(24),-1*sqrt(40),-1*sqrt(420)]
    
    gaunt_table[3][3][-2] = OffsetArray{Vector{Float64}}(undef, -2:2)
    gaunt_table[3][3][-2][-2] = [1,0,-7,6]
    gaunt_table[3][3][-2][-1] = [0,sqrt(15),sqrt(32),-1*sqrt(105)]
    gaunt_table[3][3][-2][0] = [0,-1*sqrt(20),-1*sqrt(3),4*sqrt(14)]
    gaunt_table[3][3][-2][1] = [0,0,-1*sqrt(14),-1*sqrt(378)]
    gaunt_table[3][3][-2][2] = [0,0,sqrt(70),sqrt(504)]
    
    gaunt_table[3][3][-3] = OffsetArray{Vector{Float64}}(undef, -3:3)
    gaunt_table[3][3][-3][-3] = [1,-5,3,-1]
    gaunt_table[3][3][-3][-2] = [0,5,-1*sqrt(30),sqrt(7)]
    gaunt_table[3][3][-3][-1] = [0,-1*sqrt(10),sqrt(54),-1*sqrt(28)]
    gaunt_table[3][3][-3][0] = [0,0,-1*sqrt(63),sqrt(84)]
    gaunt_table[3][3][-3][1] = [0,0,sqrt(42),-1*sqrt(210)]
    gaunt_table[3][3][-3][2] = [0,0,0,sqrt(462)]
    gaunt_table[3][3][-3][3] = [0,0,0,-1*sqrt(924)]
    
    gaunt_table
end

function gaunt_coeff_headings()
    #This function will generate a small(ish) table for the headings of the Gaunt coeff table
    # I am ssuming these will be important...
    
    gaunt_head = OffsetArray{OffsetArray{Vector{Float64}}}(undef, 0:3)
    
    # l1=s headings
    gaunt_head[0] = OffsetArray{Vector{Float64}}(undef, 0:3)
    gaunt_head[0][0] = [1,5]
    gaunt_head[0][1] = [sqrt(3),1]
    gaunt_head[0][2] = [1,sqrt(5)]
    gaunt_head[0][3] = [sqrt(7)]
    
    #l1=p headings
    gaunt_head[1] = OffsetArray{Vector{Float64}}(undef, 1:3)
    gaunt_head[1][1] = [1,5]
    gaunt_head[1][2] = [sqrt(15),sqrt(245)]
    # p-f - COMING SOON - see Condon, E.U & Shortley, G.H. (1953). The Theory of Atomic Spectra. Cambridge.
    
    #l1=d headings
    gaunt_head[2] = OffsetArray{Vector{Float64}}(undef, 2:3)
    gaunt_head[2][2] = [1,7,21]
    gaunt_head[2][3] = [sqrt(35),3*sqrt(35),6*sqrt(254)]
    
    #l1=f headings
    gaunt_head[3] = OffsetArray{Vector{Float64}}(undef, 3:3)
    gaunt_head[3][3] = [1,15,33,429/5]
    gaunt_head
end


#Crystal Field Tables - operator tables
 
function V_CFE_d_cubic()
    # This is for d electrons, it is a 5x5 matrix since only ml matters when it comes to the value of V_CEF
    # See multiplet tutorial, section 2.9 - Crystal fields for more info
    
    V_CEF = OffsetArray{Float64}(undef,-2:2,-2:2).*0
    
    # Diagonal dudes
    V_CEF[-2,-2] = 1
    V_CEF[-1,-1] = -4
    V_CEF[0,0] = 6
    V_CEF[1,1] = -4
    V_CEF[2,2] = 1
    
    # Off-diagonal dudes
    V_CEF[-2,2] = 5
    V_CEF[2,-2] = 5
    
    return V_CEF
end
