`test.cpp` contains some tests that help to ensure that modifications to the
library won't change the calculations that were done previously. At present,
the tests test for exact equality of the numbers calculated by the test and the
hard-coded values. Different compilers may apply different optimizations (or in
a different order), and due to limited precision of floating point numbers it
may result in calculated values being slightly different which will result in a
test failure, with the test complaining about inequality in a few least
significant digits. If you see that, don't worry about it &mdash; it is the
test suite that must be redesigned.

To compile the tests, Boost [test][boost.test] library is required. To compile
and run inexpensive tests, execute `make test` in the parent directory. These
tests test functions neglecting non-electromagnetic interactions. To compile
and run all tests, execute `make test_all`. Functions taking into account
non-electromagnetic interactions will take a long time to calculate and will
use all CPU cores available. The output of the tests should be
self-explanatory.

[boost.test]:  https://www.boost.org/doc/libs/1_84_0/libs/test/doc/html/index.html
