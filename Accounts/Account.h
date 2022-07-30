#ifndef ACCOUNT_H_
#define ACCOUNT_H_
#include <fprecision.h>
#include <stddef.h>
#include <stdlib.h>

namespace Accounts
{
    class A
    {
    public:
        float_precision getBalance(void);
        void setBalance(float_precision balanceToSet);

    private:
        float_precision balance;
    }
}
#endif // ACCOUNT_H_