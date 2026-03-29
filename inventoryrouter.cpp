// InventoryRouter.cpp
#include "InventoryRouter.h"
#include "GameInventory.h"
#include "ShadowInventory.h"

bool UInventoryRouter::AddItem_Implementation(TScriptInterface<IItem> Item, int32 Quantity)
{
    if (GameInventory && GameInventory->Execute_AddItem(GameInventory, Item, Quantity))
        return true;

    if (bShadowEnabled && ShadowInventory)
        return ShadowInventory->Execute_AddItem(ShadowInventory, Item, Quantity);

    return false;
}

bool UInventoryRouter::RemoveItem_Implementation(FName ItemId, int32 Quantity)
{
    if (GameInventory && GameInventory->Execute_RemoveItem(GameInventory, ItemId, Quantity))
        return true;

    if (bShadowEnabled && ShadowInventory)
        return ShadowInventory->Execute_RemoveItem(ShadowInventory, ItemId, Quantity);

    return false;
}

int32 UInventoryRouter::GetQuantity_Implementation(FName ItemId) const
{
    int32 Total = 0;
    if (GameInventory)
        Total += GameInventory->Execute_GetQuantity(GameInventory, ItemId);
    if (bShadowEnabled && ShadowInventory)
        Total += ShadowInventory->Execute_GetQuantity(ShadowInventory, ItemId);
    return Total;
}
