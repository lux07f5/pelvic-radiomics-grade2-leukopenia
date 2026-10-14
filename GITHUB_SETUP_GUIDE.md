# GitHub 注册与安装指南（写给非开发者）

> 本指南面向**没有 GitHub 账号、没用过 Git**的临床医生/科研人员。按步骤操作，**30–60 分钟**可完成全部准备。

---

## 一、GitHub 是什么？为什么学术圈要用它？

GitHub 是全球最大的**代码托管平台**（github.com），类似学术界的"论文配套代码 + 数据"标准仓库。

| 平台 | 用途 | 类比 |
|---|---|---|
| ResearchGate | 学者社交网络 | 学术朋友圈 |
| ORCID | 学者身份证 | DOI for authors |
| **GitHub** | **代码 + 数据托管** | **论文的"实验记录本"** |

**为什么 RoO 投稿要 GitHub？**
1. **可复现性（Reproducibility）**：审稿人可以直接看你的代码，验证你论文中的 AUC、P 值是不是真的
2. **加分项**：RoO 近年要求 "Code and Data Availability Statement"，有 GitHub 链接 = 审稿印象分 +20%
3. **引用率**：GitHub 仓库有 DOI 可被引用，其他研究者复用你的方法
4. **学术诚信**：公开代码让"造假"变得几乎不可能

---

## 二、注册 GitHub 账号（10 分钟）

### 步骤 1：打开 GitHub 首页
- 浏览器（Chrome / Edge / Safari）输入：**github.com**
- 点击右上角 **Sign up**（注册）

### 步骤 2：填注册信息
GitHub 注册完全免费，**不要选任何付费套餐**。

```
Username（用户名）：建议用你学术姓名的英文，如 chenweijie-2026
              ⚠️ 警告：用户名一旦确定就很难改，建议现在就想好
              （以后所有仓库 URL 都会带这个用户名）

Email address（邮箱）：建议用你的学术邮箱（如 @hospital.edu.cn）
              ⚠️ 千万不要用 QQ 邮箱或 163 邮箱（GitHub 经常屏蔽）

Password（密码）：至少 15 位，含大小写+数字+符号
              （GitHub 现在强制要求）
```

### 步骤 3：验证
- **Captcha 验证**：看图选图
- **邮箱验证**：GitHub 给你邮箱发一封邮件，点里面的链接

### 步骤 4：选套餐
- 会问你 "What kind of work are you mainly doing on GitHub?"
- 选 **"Research"** 或 **"Other"** 都可以
- **不要勾选任何付费的 Copilot / Pro 等附加项**
- 顶部会显示 **"Welcome to GitHub, free for individual users"** —— 这就对了

### 步骤 5：个性化问卷（可跳过）
- 会有 2–3 个调查问题（"How experienced are you with programming?" 选 "Beginner"）
- 最下面有个 **"Skip this step"** —— 直接跳过
- 进入个人主页（github.com/你的用户名）就成功了

---

## 三、安装 Git 客户端（10 分钟）

Git 是把本地代码推送到 GitHub 的工具。**注意：以下按平台选一种即可；鸿蒙系统另见「鸿蒙用户」小节，也可以完全跳过安装，改用 API 上传。**

### 鸿蒙用户（HarmonyOS PC / OpenHarmony）

鸿蒙有两条路，另外还有一条**完全不用装 Git** 的路：

**路 1：应用市场装 GitNext（最推荐）**
1. 打开桌面「应用市场」
2. 搜索 **GitNext** → 安装
3. 打开后配置：用户名、邮箱（须与 GitHub 注册邮箱一致）、默认仓库存储路径
4. 终端执行 `git -v`，看到 `git version 2.x.x (GitNext for HarmonyOS)` 即成功
5. 中文文件名报错的话，执行一次 `git config --global core.quotepath false`

**路 2：融合终端 `dnf` 安装（HarmonyOS 7.0+，华为官方文档写法）**

```bash
sudo dnf install git -y
git --version
```
前提是系统里能找到「融合终端 / HiShell」入口；找不到说明版本未开放。

**路 3：不装 Git，直接用 GitHub API 上传**

Git 客户端本质只是"把文件通过 HTTP 传给 GitHub"，这件事用 REST API 就能做。
本仓库配套的 `github_upload.js` 可以一次上传全部文件并生成**单个干净提交**，
只需要一个 GitHub token，不需要安装任何软件。

> ⚠️ 第三方 HarmonyBrew（`brew install git`）也可行，但属于非官方渠道且安装脚本
> 会直接执行远程代码，**有安全风险，不推荐**。

### Windows 用户

1. 打开 **https://git-scm.com/download/win**
2. 下载会自动开始（约 50 MB）
3. 双击安装，**全部按默认选项一路 Next**（不要改任何东西）
4. 安装完后，从开始菜单打开 **"Git Bash"**（一个黑窗口）

### Mac 用户

打开"终端"（Terminal.app），输入：

```bash
# 方式 1：装 Homebrew 然后装 Git（推荐）
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
brew install git

# 方式 2：装 Xcode 命令行工具（自带 Git）
xcode-select --install
```

### Linux 用户

```bash
# Ubuntu / Debian
sudo apt install git

# CentOS / RHEL
sudo yum install git
```

### 验证安装

打开终端（Windows: Git Bash；Mac: Terminal），输入：

```bash
git --version
```

应该看到类似 `git version 2.42.0` —— 成功。

---

## 四、配置 Git（5 分钟）

在终端里输入（把名字和邮箱换成你自己的）：

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@hospital.edu.cn"
```

⚠️ **邮箱必须与 GitHub 注册邮箱一致**，否则你 commit 时会报错。

---

## 五、下一步：阅读 `PUSH_GUIDE.md`

注册完 GitHub 账号 + 安装好 Git 后，看仓库根目录的 **`PUSH_GUIDE.md`**，按里面的 5 步把代码推上去。

---

## 常见问题

### Q1: 我能不能用学校的 VPN / 校园网？
可以。如果连不上 GitHub，看 `Q5`。

### Q2: 我能不能用 QQ 邮箱 / 163 邮箱？
**不推荐**。很多 GitHub 通知邮件会被国内邮箱误判为垃圾邮件。如果只有国内邮箱，去申请一个学术邮箱（如 Gmail、Outlook、机构邮箱）。

### Q3: 用户名能不能用中文？
**不能**。用户名只支持英文、数字、连字符。建议 `姓拼音+名拼音首字母` 格式，如 `chenwj2026`。

### Q4: Git 是什么？和 GitHub 是什么关系？
- **Git**：本地装的一款软件，用来"管理代码版本"（类似文档历史记录）
- **GitHub**：网上的一个网站，用来"存放 Git 管理的代码"

**Git = 你的本地仓库，GitHub = 云端备份 + 分享**。

### Q5: 连不上 github.com 怎么办？
- 检查网络：可以打开网页吗？
- 如果网页能打开但 `git push` 很慢：尝试开启 VPN（GitHub 服务器在国外）
- 如果是公司 / 医院内网：可能禁用了 22 端口，联系网管开端口或用 HTTPS（见 PUSH_GUIDE）

### Q6: 我完全不想注册 GitHub，能直接投稿吗？
**能**，但在 "Data Availability Statement" 写：
> "The analysis code is available from the corresponding author upon reasonable request."

缺点是审稿人可能懒得联系你 → 论文评分降低。**强烈建议建 GitHub**。

---

## 遇到问题？

1. 看错误信息：GitHub/Git 的错误信息通常很明确（英文）
2. 复制错误信息前 100 字 + 你的操作步骤，发给 AI 助手
3. 不要慌：90% 的问题都在这 2 个文档里

---

*预计完成时间：30–60 分钟。下一步：`PUSH_GUIDE.md`*
