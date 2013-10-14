//NumericAlphabet.h
//
//Alphabets whose elements are integers.

#ifndef NumericAlphabet_h
#define NumericAlphabet_h

class NumericAlphabet
{
    private:
        int Length;

    public:
        NumericAlphabet(int n = 0) { assert(n >= 0); set_length(n); }
        ~NumericAlphabet() {};

        void set_length(int n) { Length = n; }
        int length(void) const { return Length; }

        int index(int i) const { assert(0 <= i  &&  i < Length); return i; }
        int undindex(int i) const { assert(0 <=i && i < Length); return i; }
};

#endif
