#ifndef ACCOUNT_H_
#define ACCOUNT_H_
#include <fprecision.h>
#include <stddef.h>
#include <stdlib.h>
#include "expense.h"

namespace Accounts
{
    class Account
    {
    public:
        float_precision getBalance(void);
        void setBalance(float_precision balanceToSet);
        virtual void deductExpense(Expenses::Expense expense) = 0;

    private:
        float_precision balance;
    };
}
#endif // ACCOUNT_H_