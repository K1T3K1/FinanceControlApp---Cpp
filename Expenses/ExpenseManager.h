#ifndef EXPENSEMANAGER_H_
#define EXPENSEMANAGER_H_

#include "Expense.h"
#include "Account.h"

namespace Expenses
{
    class ExpenseManager
    {
    public:
        static void addExpense(Accounts::Account *recAccount, std::string expenseTitle, float_precision expenseCost, tm *expenseDate);

    private:
    };
}

#endif // EXPENSEMANAGER_H_