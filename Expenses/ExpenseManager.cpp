#include "Expense.h"
#include "ExpenseManager.h"

namespace Expenses
{
    void ExpenseManager::addExpense
    (Accounts::Account *recAccount, std::string expenseTitle, float_precision expenseCost, tm *expenseDate)
    {
        Expense* ExpensetoAdd = new Expense(expenseTitle, expenseCost, expenseDate);
        *recAccount->expenseToAdd(ExpensetoAdd);

    }
}