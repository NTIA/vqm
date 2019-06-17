function offset = find_test(test_structs, test)
% FIND_ORIGINAL
%  Find a test structure's offset, given the test name
% SYNTAX
%  offset = find_test(test_structs, test);
% DESCRIPTION
%  When offset = find_test(test_structs, test); is run,
%  'test_structs' is an array of tests, all of the same format as GTests.
%
%  'test' is the name of the test you want to find.  Function 'find_test' 
%  locates and returns the array address of the requested test.  If the 
%  test cannot be found, and error will occur.

all_test = [test_structs.name];
[temp, offset] = find(strcmp(test,all_test));

if length(offset) == 0,
    error('Test %s cannot be found', test);
end
