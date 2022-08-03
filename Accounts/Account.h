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
        float_precision getBalance(void);
        void setBalance(float_precision balanceToSet);
        void expenseToAdd(Expenses::Expense expense);
        virtual void deductExpense(Expenses::Expense *expense) = 0;

    private:
        float_precision balance;
        std::vector<Expenses::Expense *> associatedExpenses;
    };
}
#endif // ACCOUNT_H_