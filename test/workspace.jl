
function shift(vec, k)
temp = vec[1 : end - k]
vec[1 : k] = 0
vec[1 + k : end] = temp
end

a = [1; 1; [0 for j = 1 : 8]]

for k = 1 : 2
    a = shift(a, 1)
    println(a)
end
