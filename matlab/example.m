struct1=struct('a',zeros(10,1),'b',zeros(100,1),'c',uint32(0));
struct_array1(1:100,1)=struct1;
parfor i=1:100
    a=rand(10,1);
    b=rand(100,1);
    c=uint32((a(1)+b(1))*1000);
    struct_array1(i).a=a;
    struct_array1(i).b=b;
    struct_array1(i).c=c;
end