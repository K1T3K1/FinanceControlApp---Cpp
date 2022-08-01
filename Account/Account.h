#ifndef ACCOUNT_H_
#define ACCOUNT_H_
#include <fprecision.h>
#include <stddef.h>
#include <stdlib.h>
#include "Expense.h"

namespace Accounts
{
    class Account
    {
    public:
        Account(std::string Name, float_precision startBalance);
        float_precision getBalance(void);
        void setBalance(float_precision balanceToSet);
        virtual void deductExpense(Expenses::Expense expense) = 0;

    private:
        std::string accountSignature;
        float_precision balance;
        std::vector<Expenses::Expense*> latestExpenses;
    };
}
#endif // ACCOUNT_H_