#ifndef ACCOUNTMANAGER_H_
#define ACCOUNTMANAGER_H_
#include "Account.h"
#include "Expense.h"

namespace Accounts
{
    class AccountManager
    {
    public:
        void addAccount();
        Account getAccount(int index);

    private:
        std::vector<Account *> AccountList;
    };
}

#endif // ACCOUNTMANAGER_H_