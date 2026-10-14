# 推送指南：5 步把代码推上 GitHub

> **前提**：你已经按 `GITHUB_SETUP_GUIDE.md` 注册了账号 + 安装了 Git。
> 本指南：**Windows 用户用 Git Bash，Mac 用户用 Terminal**。命令几乎一样（路径分隔符不同）。

---

## 步骤 0：在 GitHub 网站上创建空仓库（3 分钟）

1. 浏览器打开 **github.com**，登录
2. 右上角 **+** → **New repository**
3. 填表：

| 字段 | 填什么 |
|---|---|
| **Repository name** | `pelvic-radiomics-grade2-leukopenia` |
| **Description** | `Planning CT radiomics model for Grade ≥ 2 leukopenia after pelvic radiotherapy (60 patients, 15 events)` |
| **Public / Private** | **Public**（公开） |
| **Initialize this repository with** | **全部不勾选**（README/.gitignore/LICENSE 都不勾——我们本地已经有了） |

4. 点 **Create repository**
5. 创建后页面会显示一段命令提示，**先别关这个页面**，下一步要用里面的 URL

---

## 步骤 1：把仓库文件夹搬到你的电脑上（5 分钟）

把整个 `github_repo/` 文件夹复制到你电脑上任意位置，建议放桌面或文档目录：

```
C:\Users\你的用户名\Desktop\pelvic-radiomics-grade2-leukopenia\      (Windows)
/Users/你的用户名/Desktop/pelvic-radiomics-grade2-leukopenia/        (Mac)
```

⚠️ **不要放在有中文路径的目录里**（如 `桌面\我的论文\`），Git 对中文路径偶尔有问题。

---

## 步骤 2：打开终端，进入仓库目录

### Windows（Git Bash）

右键点击 `pelvic-radiomics-grade2-leukopenia` 文件夹 → **Git Bash Here**

或者：
- 打开 Git Bash
- 输入 `cd /c/Users/你的用户名/Desktop/pelvic-radiomics-grade2-leukopenia`

### Mac（Terminal）

打开 Terminal.app，输入：
```bash
cd /Users/你的用户名/Desktop/pelvic-radiomics-grade2-leukopenia
```

验证：输入 `ls`（Mac/Linux）或 `dir`（Windows），应能看到 `README.md`、`LICENSE`、`.gitignore` 等文件。

---

## 步骤 3：5 条 Git 命令（一行一行复制粘贴执行）

```bash
# 1. 初始化 Git 仓库（只在第一次推送时需要）
git init

# 2. 把主分支改名为 main（GitHub 默认是 main，老版本是 master）
git branch -M main

# 3. 把所有文件加入暂存区
git add .

# 4. 首次提交（必须写提交信息）
git commit -m "Initial release: 60-patient pure RT radiomics prediction model"

# 5. 关联到你的 GitHub 仓库
#    把下面 URL 里的 [你的用户名] 替换成你 GitHub 注册时的用户名
git remote add origin https://github.com/[你的用户名]/pelvic-radiomics-grade2-leukopenia.git

# 6. 推送到 GitHub（这一步会要求你输入 GitHub 用户名 + 密码/token）
git push -u origin main
```

⚠️ **第 6 步可能要求输入密码**：
- 2021 年起，GitHub 不再支持直接用密码推送
- 你需要创建一个 **Personal Access Token (PAT)**
- 创建方法：github.com → 右上角头像 → Settings → Developer settings → Personal access tokens → Tokens (classic) → Generate new token
- 选 **`repo`** 权限，过期选 90 天或 No expiration
- 把生成的 token 当作密码粘贴（token 通常以 `ghp_` 开头）

---

## 步骤 4：验证推送成功

浏览器打开：
```
https://github.com/[你的用户名]/pelvic-radiomics-grade2-leukopenia
```

应该能看到：
- README.md 渲染成漂亮的网页
- LICENSE 显示 MIT 徽章
- 文件树显示 scripts/、data/、results/ 三个目录

如果看到了，**恭喜！推送成功** ✅

---

## 步骤 5：把 URL 填回论文

打开 `论文_Methods_EN.md`（你的论文英文版），找到 §2.10 "Statistical software" 段，把：

```
GitHub 仓库链接 [待补]
```

替换为：

```
https://github.com/[你的用户名]/pelvic-radiomics-grade2-leukopenia
```

同时在论文末尾的 "Data Availability Statement" 段写：

> Data Availability: The de-identified analysis code and summary statistics are publicly available at https://github.com/[你的用户名]/pelvic-radiomics-grade2-leukopenia. Raw imaging data are not publicly available due to patient privacy concerns but may be made available to qualified researchers upon reasonable request to the corresponding author, subject to IRB approval.

---

## 常见错误速查

| 错误信息 | 原因 | 解决 |
|---|---|---|
| `fatal: not a git repository` | 当前目录不是 git 仓库 | 确认 `cd` 到了正确的文件夹；或重新 `git init` |
| `Permission denied (publickey)` | 用了 SSH 但没配公钥 | 改用 HTTPS：`git remote set-url origin https://...` |
| `Authentication failed` | 密码错了 / token 失效 | 重新生成 PAT（步骤 3 末尾） |
| `Updates were rejected because...` | GitHub 仓库已有内容（README等） | 第一次用我们这个本地仓库推，**确保步骤 0 时不勾选任何初始化选项** |
| `filename too long` / 中文文件名错误 | Windows + 中文路径 | 把仓库移到纯英文路径下 |
| 推送后网页显示 404 | URL 拼错了 | 检查 `[你的用户名]` 是否替换正确 |
| 推送后网页没显示某些文件 | 文件被 `.gitignore` 过滤了 | `cat .gitignore` 看规则；调整后 `git add . -f` 强制添加（⚠️ 不要加 patient 文件） |

---

## 后续操作：怎么更新代码？

你以后每次修改代码后，只需要 3 条命令：

```bash
git add .                                # 把改过的文件加入暂存区
git commit -m "简短描述你改了什么"        # 提交
git push                                  # 推送到 GitHub
```

例如：
```bash
git add scripts/python/make_figures.py
git commit -m "Update Figure 5 with corrected axis labels"
git push
```

---

## 怎么删除不小心推上去的隐私文件？

**紧急操作**（如果误推了含患者姓名的文件）：

```bash
# 1. 从仓库历史中删除
git rm --cached <文件名>
git commit -m "Remove accidentally uploaded PHI file"
git push

# 2. 同步删除 GitHub 历史（用 BFG Repo-Cleaner，专业级工具）
#    见：https://rtyley.github.io/bfg-repo-cleaner/
#    警告：操作前备份本地仓库
```

**重要**：如果真的误推了 PHI，立刻联系 GitHub Support 申请删除（https://support.github.com/contact）—— Git 仓库历史即使删除文件，commit 记录里仍可恢复。

---

## 想切换为 SSH 而非 HTTPS？（可选，进阶）

HTTPS 每次推送要输 token（或配 credential helper）。SSH 配置好后免密码。

详细见：https://docs.github.com/en/authentication/connecting-to-github-with-ssh

简单版：
```bash
# 生成 SSH key（默认位置 ~/.ssh/id_ed25519）
ssh-keygen -t ed25519 -C "your.email@hospital.edu.cn"
# 一路回车，不设密码

# 把公钥复制到剪贴板
# Mac:
pbcopy < ~/.ssh/id_ed25519.pub
# Windows (Git Bash):
clip < ~/.ssh/id_ed25519.pub

# 然后到 GitHub 网页：Settings → SSH and GPG keys → New SSH key → 粘贴

# 改 remote URL 为 SSH
git remote set-url origin git@github.com:[你的用户名]/pelvic-radiomics-grade2-leukopenia.git

# 测试
ssh -T git@github.com
```

---

## 完成检查清单

- [ ] github.com 上能看到仓库页面
- [ ] README.md 渲染正常
- [ ] 文件树显示完整
- [ ] 论文 Methods §2.10 + Data Availability 段已填入 URL
- [ ] 关闭 .git 之前没有警告（说明没有 PHI 被过滤掉）

---

*预计完成时间：20–30 分钟（不算下载大文件）*
