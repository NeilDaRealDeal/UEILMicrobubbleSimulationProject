% a = 0
% for i = 1:1e11
%     a = a+1;
% end
% disp(a)
a = 0
parfor i = 1:1e11
    a = a+1;
end
disp(a)

