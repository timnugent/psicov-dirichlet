#ifndef IOSMACROS_H
#define IOSMACROS_H

// These macros provide the left, right, ... manipulators even
// for old compilers that lack the new iomanips.

#include <iomanip>
#define IOSLEFT resetiosflags(ios::adjustfield) << setiosflags(ios::left)
#define IOSRIGHT resetiosflags(ios::adjustfield) << setiosflags(ios::right)
#define IOSFIXED resetiosflags(ios::floatfield) << setiosflags(ios::fixed)
#define IOSSCIENTIFIC resetiosflags(ios::floatfield) << setiosflags(ios::scientific)

using namespace std;

#endif
