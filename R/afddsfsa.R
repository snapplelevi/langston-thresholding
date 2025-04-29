# 
# inlist <- c(10,20,30,40,50,50,60,60,70,70,80,80,90,10,11,11,12,12)
# 
# 
# l <- (function(input_list){
#     
#     dup_list <- duplicated(input_list)
#     
#     # First element of list will always be unique, so start at the second
#     # for the while loop logic to be clean.
#     rep_list <- c(1)
#     ind <- 2
#     inner_ind <- 2
#     while(inner_ind <= length(dup_list)){
#         
#         streak <- 1
#         
#         while(inner_ind < length(dup_list) && dup_list[[inner_ind+1]] == TRUE){
#             streak <- streak + 1
#             inner_ind <- inner_ind + 1
#         }
#         
#         rep_list <- append(rep_list, rep(ind, streak))
#         ind <- ind + 1
#         inner_ind <- inner_ind + 1
#     }
#     rep_list
#     
# })(inlist)
# 
# 
# unique_list <- unique(inlist) * 2
# 
# print(inlist)
# print(unique_list)
# 
# print(l)
# 
# 
# print("comparing:::")
# print("original first, then replist thing")
# print(inlist)
# print(length(inlist))
# repped <- unique_list[l]
# print(repped)
# print(length(repped))
